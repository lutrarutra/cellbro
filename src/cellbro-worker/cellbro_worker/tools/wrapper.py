import json
from typing import Callable, Sequence

from cellbro_db import types


from .. import celery_app
from . import worker_redis
from .OutputCaptureHandler import StdoutCaptureHandler

def worker_task(
    task_name: str,
    read_resources: Sequence[str],
    write_resources: Sequence[str] | None = None,
    complete_steps: types.ChecklistStep | list[types.ChecklistStep] | None = None,
    trigger_events: str | list[str] | None = None,
    notify: bool = False,
) -> Callable:
    if complete_steps is None:
        complete_steps = []
    elif isinstance(complete_steps, types.ChecklistStep):
        complete_steps = [complete_steps]

    if trigger_events is None:
        trigger_events = []
    elif isinstance(trigger_events, str):
        trigger_events = [trigger_events]

    def decorator(func: Callable) -> Callable:
        @celery_app.task(name=task_name, bind=True)
        def wrapper(*args, **kwargs):
            r_res = set(read_resources)
            w_res = set(write_resources) if write_resources else set()
            all_requested = r_res | w_res

            with celery_app.state_lock:
                while (
                    any(res in celery_app.active_writers for res in all_requested) or
                    any(res in celery_app.active_readers for res in w_res)
                ):
                    celery_app.state_lock.wait()

                if w_res:
                    celery_app.active_writers.update(w_res)
                celery_app.active_readers.update(r_res)

            worker_redis.set("current_task", task_name)
            worker_redis.publish("current_task", task_name)
            worker_redis.publish("task_started", task_name)
            if notify:
                worker_redis.publish("notify", json.dumps({"message": f"Task {task_name} started.", "category": "info"}))

            try:
                with StdoutCaptureHandler("general", worker_redis):
                    result = func(*args, **kwargs)
                    print(f"Task {task_name} completed.")
                    
                worker_redis.publish("task_completed", task_name)

                for complete_step in complete_steps:
                    worker_redis.set(f"step:{complete_step}", "completed")
                    worker_redis.publish("step_completed", complete_step)

                for event in trigger_events:
                    worker_redis.publish("event_triggered", event)

                if notify:
                    worker_redis.publish("notify", json.dumps({"message": f"Task {task_name} completed.", "category": "success"}))

            except Exception as e:
                print(f"An error occurred in task {task_name}: {e}")
                result = None
            finally:
                with celery_app.state_lock:
                    if w_res:
                        celery_app.active_writers.difference_update(w_res)
                    celery_app.active_readers.difference_update(r_res)
                    celery_app.state_lock.notify_all()

            worker_redis.set("current_task", "idle")
            worker_redis.publish("current_task", "idle")
            return result
        return wrapper
    return decorator