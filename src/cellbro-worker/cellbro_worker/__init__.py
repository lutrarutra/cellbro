import os
import threading
from celery import Celery, Task
from celery.signals import worker_ready
import anndata as ad
from redis import Redis

from cellbro_db import DBHandler, types


REDIS_PORT = int(os.environ.get("REDIS_PORT", 6379))

class CellBroWorker(Celery):
    def __init__(self, name: str, redis_port: int = 6379):
        super().__init__(name, broker=f"redis://redis-cache:{redis_port}/4")
        self.adata: ad.AnnData = None  # type: ignore
        self.state_lock = threading.Condition()
        self.active_readers = set()
        self.active_writers = set()
        self.running_task_count = 0
        self.r = Redis(host="redis-cache", port=int(os.environ["REDIS_PORT"]), db=5, decode_responses=True)

    
class CellBroTask(Task):
    @property
    def app(self) -> CellBroWorker:
        return super().app
    
    @property
    def adata(self) -> ad.AnnData:
        return self.app.adata
    
    @adata.setter
    def adata(self, value: ad.AnnData):
        self.app.adata = value

    @property
    def r(self) -> Redis:
        return self.app.r
    
    @property
    def db(self) -> DBHandler:
        db = DBHandler(auto_commit=True)
        db.connect(
            user=os.environ["POSTGRES_USER"],
            password=os.environ["POSTGRES_PASSWORD"],
            host="postgres",
            port=os.environ["POSTGRES_PORT"],
            db=os.environ["POSTGRES_DB"],
        )
        return db
    
    @property
    def task_id(self) -> str:
        return self.request.id



celery_app = CellBroWorker("cellbro-worker", redis_port=REDIS_PORT)
celery_app.conf.update(
    task_track_started=True,
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone=os.environ.get("TZ", "UTC"),
)

@worker_ready.connect
def cleanup_on_start(sender, **kwargs):
    db = DBHandler(auto_commit=True)
    db.connect(
        user=os.environ["POSTGRES_USER"],
        password=os.environ["POSTGRES_PASSWORD"],
        host="postgres",
        port=os.environ["POSTGRES_PORT"],
        db=os.environ["POSTGRES_DB"],
    )
    for step in types.ChecklistStep:
        celery_app.r.delete(f"step:{step}")


from . import tasks
celery_app.autodiscover_tasks()