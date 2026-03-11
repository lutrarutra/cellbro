import os
from taskiq import (
    TaskiqEvents, 
    TaskiqState, 
    TaskiqMiddleware, 
    TaskiqMessage, 
    TaskiqResult, 
)
from taskiq_redis import RedisStreamBroker, RedisAsyncResultBackend
from cellbro_db import DBHandler

from .tools.MessageRedis import SyncMessageRedis

from . import REDIS_PORT, REDIS_URL, task_context, TaskContext

class CellBroWorkerState(TaskiqState):
    db: DBHandler

class TaskLifeCycleMiddleware(TaskiqMiddleware):
    def __init__(self):
        self.message_redis = SyncMessageRedis()
        self.message_redis.connect(host="redis-cache", port=REDIS_PORT, db=2)
        
    async def pre_execute(self, message: TaskiqMessage) -> TaskiqMessage:
        task_context.set(TaskContext(
            task_id=message.task_id,
            task_name=message.task_name,
            message_client=self.message_redis
        ))
        self.message_redis.log(message.task_id, f"Task {message.task_name} started.", category="log")
        return message
    
    async def post_execute(self, message: TaskiqMessage, result: TaskiqResult) -> None:
        if result.is_err:
            self.message_redis.error(message.task_id, f"Task {message.task_name} failed with error: {result.unwrap_err()}")
        else:
            self.message_redis.log(message.task_id, f"Task {message.task_name} completed successfully.", category="log")
        self.message_redis.complete_task(message.task_id)
    
result_backend = RedisAsyncResultBackend(redis_url=REDIS_URL.format(port=REDIS_PORT, db=4))
task_broker = RedisStreamBroker(
    url=REDIS_URL.format(port=REDIS_PORT, db=3)
).with_result_backend(result_backend).with_middlewares(TaskLifeCycleMiddleware())

@task_broker.on_event(TaskiqEvents.WORKER_STARTUP)  # type: ignore
async def startup(state: CellBroWorkerState):
    state.db = DBHandler(auto_commit=True)
    print("Connecting to database...")
    state.db.connect(
        user=os.environ["POSTGRES_USER"],
        password=os.environ["POSTGRES_PASSWORD"],
        host="postgres",
        port=os.environ["POSTGRES_PORT"],
        db=os.environ["POSTGRES_DB"],
    )

@task_broker.on_event(TaskiqEvents.WORKER_SHUTDOWN)  # type: ignore
async def cleanup_on_shutdown(state: CellBroWorkerState):
    state.db.close_connection()
    del state.db

from .tasks import io, qc