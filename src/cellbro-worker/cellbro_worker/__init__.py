import os
from celery import Celery

import anndata as ad

REDIS_PORT = int(os.environ.get("REDIS_PORT", 6379))

class CellBroWorker(Celery):
    def __init__(self, name: str, redis_port: int = 6379):
        super().__init__(name, broker=f"redis://redis-cache:{redis_port}/4")
        self.adata: ad.AnnData

celery_app = CellBroWorker("cellbro-worker", redis_port=REDIS_PORT)
celery_app.conf.update(
    task_track_started=True,
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone=os.environ.get("TZ", "UTC"),
)

from . import tasks
celery_app.autodiscover_tasks()