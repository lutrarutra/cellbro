import os
from redis.asyncio import Redis as AsyncRedis
from cellbro_db import DBHandler, SyncSession
from taskiq import TaskiqEvents, TaskiqState, TaskiqMiddleware, TaskiqMessage, TaskiqResult
from taskiq_redis import RedisStreamBroker, RedisAsyncResultBackend


REDIS_PORT = int(os.environ.get("REDIS_PORT", 6379))
result_backend = RedisAsyncResultBackend(
    redis_url=f"redis://redis-cache:{REDIS_PORT}/4",
)

class CellBroWorkerState(TaskiqState):
    db: DBHandler

class TaskLifeCycleMiddleware(TaskiqMiddleware):
    def __init__(self, redis_url: str):
        self.redis: AsyncRedis = AsyncRedis.from_url(redis_url)

    async def pre_execute(self, message: TaskiqMessage) -> TaskiqMessage:
        channel = f"task:{message.task_id}"
        await self.redis.publish(channel, "started")
        return message
    
    async def post_execute(self, message: TaskiqMessage, result: TaskiqResult) -> None:
        channel = f"task_events:{message.task_id}"
        
        if result.is_err:
            await self.redis.publish(
                channel, 
                f"Task {message.task_id} completed unsuccessfully. Error: {result.error}"
            )
        else:
            await self.redis.publish(
                channel, 
                f"Task {message.task_id} completed successfully."
            )
    
broker = RedisStreamBroker(
    url=f"redis://redis-cache:{REDIS_PORT}/4",
).with_result_backend(result_backend).with_middlewares(TaskLifeCycleMiddleware(redis_url=f"redis://redis:{REDIS_PORT}/4"))


@broker.on_event(TaskiqEvents.WORKER_STARTUP)  # type: ignore
async def startup(state: CellBroWorkerState):
    state.db = DBHandler(auto_commit=True)
    state.db.connect(
        user=os.environ["POSTGRES_USER"],
        password=os.environ["POSTGRES_PASSWORD"],
        host="postgres",
        port=os.environ["POSTGRES_PORT"],
        db=os.environ["POSTGRES_DB"],
    )

@broker.on_event(TaskiqEvents.WORKER_SHUTDOWN)  # type: ignore
async def cleanup_on_shutdown(state: CellBroWorkerState):
    state.db.close_connection()
    del state.db


def db_session(state: CellBroWorkerState) -> SyncSession:
    return state.db.open_session()


from .tasks import io, qc

    