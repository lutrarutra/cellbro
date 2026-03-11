import contextlib
import functools
from typing import Callable

from .. import current_task_id, redis_client
from .OutputCaptureHandler import StdoutCaptureHandler
from .WarningCaptureHandler import WarningCaptureHandler
from cellbro_db import types

running_task_count = 0

def step_wrapper(
    complete_steps: types.ChecklistStep | list[types.ChecklistStep] | None = None,
    trigger_events: str | list[str] | None = None,
    capture_stdout: bool = True, capture_warnings: bool = True,
    notify: bool = True
) -> Callable:
    
    complete_steps = [complete_steps] if isinstance(complete_steps, types.ChecklistStep) else (complete_steps or [])
    trigger_events = [trigger_events] if isinstance(trigger_events, str) else (trigger_events or [])
    
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            result = None
            try:
                stack = contextlib.ExitStack()
                task_id = current_task_id.get()
                r = redis_client.get()
                
                if capture_stdout:
                    stack.enter_context(StdoutCaptureHandler(task_id, r))
                    
                if capture_warnings:
                    stack.enter_context(WarningCaptureHandler(task_id, r))

                result = await func(*args, **kwargs)
                if notify:
                    r.notify(f"Task {task_id} completed")
                    
                for complete_step in complete_steps:
                    r.step(complete_step, "completed")
                    
                for event in trigger_events:
                    r.event(event)
                    
            except Exception as e:
                print(f"An error occurred in function {func.__name__}: {e}", flush=True)
                raise e
                
            return result
        return wrapper
    return decorator