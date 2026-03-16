import contextlib
import functools
from typing import Callable

from cellbro_db import types

from .. import task_context
from .OutputCaptureHandler import StdoutCaptureHandler
from .WarningCaptureHandler import WarningCaptureHandler

def worker_task(
    complete_steps: types.ChecklistStep | list[types.ChecklistStep] | None = None,
    trigger_events: str | list[str] | None = None,
    capture_stdout: bool = True, 
    capture_warnings: bool = True,
    notify: bool = True,
) -> Callable:
    
    complete_steps = [complete_steps] if isinstance(complete_steps, types.ChecklistStep) else (complete_steps or [])
    trigger_events = [trigger_events] if isinstance(trigger_events, str) else (trigger_events or [])
    
    def decorator(func: Callable) -> Callable:
        @functools.wraps(func)
        async def wrapper(*args, **kwargs):
            result = None
            _task_context = task_context.get()
            task_name = _task_context.task_name
            task_id = _task_context.task_id
            r = _task_context.message_client

            try:
                stack = contextlib.ExitStack()
                
                if capture_stdout:
                    stack.enter_context(StdoutCaptureHandler(task_id, r))
                    
                if capture_warnings:
                    stack.enter_context(WarningCaptureHandler(task_id, r))
                result = await func(*args, **kwargs)
                if notify:
                    r.notify(f"Task {task_name} Completed!", category="success")
                    
                for complete_step in complete_steps:
                    r.step(complete_step, "completed")
                    
                for event in trigger_events:
                    r.event(event)
                    
            except Exception as e:
                if notify:
                    r.notify(f"Task {task_name} Failed..", category="error")
                print(f"An error occurred in function {func.__name__}: {e}", flush=True)
                raise e
                
            return result
        return wrapper
    return decorator