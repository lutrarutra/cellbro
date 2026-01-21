import os
from celery import Celery

REDIS_PORT = int(os.environ["REDIS_PORT"])

celery_app = Celery("cellbro-tasks", broker=f"redis://redis-cache:{REDIS_PORT}/4",)

celery_app.conf.update(
    task_track_started=True,
    task_serializer="json",
    accept_content=["json"],
    result_serializer="json",
    timezone=os.environ.get("TZ", "UTC"),
    task_default_priority=5,
    task_queue_max_priority=10,
)

from . import tasks
celery_app.autodiscover_tasks()