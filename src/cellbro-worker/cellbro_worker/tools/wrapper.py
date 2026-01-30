import json
from typing import Callable, Sequence

from cellbro_db import types


from .. import celery_app, CellBroTask
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
        @celery_app.task(name=task_name, bind=True, base=CellBroTask)
        def wrapper(self: CellBroTask, *args, **kwargs):
            celery_app.running_task_count += 1
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

            self.r.set("running_tasks", celery_app.running_task_count)
            self.r.publish(self.task_id, "started")
            self.r.publish("task_started", task_name)

            if notify:
                self.r.publish("notify", json.dumps({"message": f"Task {task_name} started.", "category": "info"}))

            try:
                with StdoutCaptureHandler("general", self.r):
                    result = func(self, *args, **kwargs)
                    print(f"Task {task_name} completed.")
                    
                self.r.publish("task_completed", task_name)
                self.r.publish(self.task_id, "completed")

                for complete_step in complete_steps:
                    self.r.set(f"step:{complete_step}", "completed")
                    self.r.publish("step_completed", complete_step)

                for event in trigger_events:
                    self.r.publish("event_triggered", event)

                if notify:
                    self.r.publish("notify", json.dumps({"message": f"Task {task_name} completed.", "category": "success"}))

            except Exception as e:
                print(f"An error occurred in task {task_name}: {e}")
                self.r.publish("task_completed", task_name)
                self.r.publish(self.task_id, "failed")

                if notify:
                    self.r.publish("notify", json.dumps({"message": f"Task {task_name} failed: {e}", "category": "error"}))

                result = None
            finally:
                with celery_app.state_lock:
                    if w_res:
                        celery_app.active_writers.difference_update(w_res)
                    celery_app.active_readers.difference_update(r_res)
                    celery_app.state_lock.notify_all()

                celery_app.running_task_count -= 1
                self.r.set("running_tasks", celery_app.running_task_count)

            return result
        return wrapper
    return decorator