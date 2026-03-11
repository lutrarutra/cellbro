import os
from typing import TYPE_CHECKING
from contextvars import ContextVar

from taskiq import (
    TaskiqDepends,
    Context as TaskiqContext,
)

from cellbro_db import SyncSession

if TYPE_CHECKING:
    from .tools.MessageRedis import SyncMessageRedis

REDIS_PORT = int(os.environ.get("REDIS_PORT", 6379))
REDIS_URL = "redis://redis-cache:{port}/{db}"


class TaskContext:
    task_id: str
    task_name: str
    message_client: "SyncMessageRedis"

    def __init__(self, task_id: str, task_name: str, message_client: "SyncMessageRedis"):
        self.task_id = task_id
        self.task_name = task_name
        self.message_client = message_client

task_context: ContextVar[TaskContext] = ContextVar("task_context", default=TaskContext(task_id="", task_name="", message_client=None))  # type: ignore

def db_session(context: TaskiqContext = TaskiqDepends()) -> SyncSession:
    state: CellBroWorkerState = context.state  # type: ignore
    return state.db.open_session()