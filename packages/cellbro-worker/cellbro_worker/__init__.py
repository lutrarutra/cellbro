import os
from typing import TYPE_CHECKING
from contextvars import ContextVar

from taskiq import (
TaskiqDepends,
    Context
)

from cellbro_db import SyncSession

if TYPE_CHECKING:
    from .tools.MessageRedis import SyncMessageRedis

REDIS_PORT = int(os.environ.get("REDIS_PORT", 6379))
REDIS_URL = "redis://redis-cache:{port}/{db}"

current_task_id: ContextVar[str] = ContextVar("current_task_id", default="local_run")
redis_client: ContextVar["SyncMessageRedis"] = ContextVar("redis_client", default=None)  # type: ignore

def db_session(context: Context = TaskiqDepends()) -> SyncSession:
    state: CellBroWorkerState = context.state  # type: ignore
    return state.db.open_session()