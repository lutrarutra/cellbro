from typing import Callable

from cellbro_db import types

from .. import celery_app
from . import worker_redis
from .OutputCaptureHandler import StdoutCaptureHandler

def worker_task(
    task_name: str,
    complete_steps: types.ChecklistStep | list[types.ChecklistStep] | None = None,
    trigger_events: str | list[str] | None = None,
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
            worker_redis.set("current_task", task_name)
            worker_redis.publish("current_task", task_name)
            worker_redis.publish("task_started", task_name)
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

            except Exception as e:
                print(f"An error occurred in task {task_name}: {e}")
                result = None
            worker_redis.set("current_task", "idle")
            worker_redis.publish("current_task", "idle")
            return result
        return wrapper
    return decorator