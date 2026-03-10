import json
import functools
import redis
from typing import Callable, Sequence

# Import your Taskiq broker
from .. import broker
from .OutputCaptureHandler import StdoutCaptureHandler
from cellbro_db import types

# Assuming you have a centralized Redis connection URL
REDIS_URL = "redis://redis-cache:6377/4"

running_task_count = 0

def worker_task(
    task_name: str,
    read_resources: Sequence[str],
    write_resources: Sequence[str] | None = None,
    complete_steps: types.ChecklistStep | list[types.ChecklistStep] | None = None,
    trigger_events: str | list[str] | None = None,
    notify: bool = False,
) -> Callable:
    
    # Standardize inputs
    complete_steps = [complete_steps] if isinstance(complete_steps, types.ChecklistStep) else (complete_steps or [])
    trigger_events = [trigger_events] if isinstance(trigger_events, str) else (trigger_events or [])
    
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        def wrapper(*args, **kwargs):
            global running_task_count
            
            r = redis.from_url(REDIS_URL)
            
            task_id = task_name 
            
            r_res = set(read_resources)
            w_res = set(write_resources) if write_resources else set()
            all_requested = r_res | w_res
                
            running_task_count += 1
            r.set("running_tasks", running_task_count)
            
            r.publish(task_id, "started")
            r.publish("task_started", task_name)
            if notify:
                r.publish("notify", json.dumps({"message": f"Task {task_name} started.", "category": "info"}))
                
            result = None
            try:
                with StdoutCaptureHandler("general", r):
                    result = func(*args, **kwargs)
                    print(f"Task {task_name} completed.")
                    
                r.publish("task_completed", task_name)
                r.publish(task_id, "completed")
                
                for complete_step in complete_steps:
                    r.set(f"step:{complete_step}", "completed")
                    r.publish("step_completed", complete_step)
                    
                for event in trigger_events:
                    r.publish("event_triggered", event)
                    
                if notify:
                    r.publish("notify", json.dumps({"message": f"Task {task_name} completed.", "category": "success"}))
                    
            except Exception as e:
                print(f"An error occurred in task {task_name}: {e}")
                r.publish("task_completed", task_name)
                r.publish(task_id, "failed")
                if notify:
                    r.publish("notify", json.dumps({"message": f"Task {task_name} failed: {e}", "category": "error"}))
                raise e # ⚠️ Re-raise so Taskiq knows it failed!
                
            finally:
                running_task_count -= 1
                r.set("running_tasks", running_task_count)
                
            return result
        return broker.task(task_name=task_name)(wrapper)

    return decorator