import os
from cellbro_db import DBHandler
from redis import Redis

from .. import celery_app
from ..tools.OutputCaptureHandler import StdoutCaptureHandler

stdout_redis = Redis(host="redis-cache", port=int(os.environ["REDIS_PORT"]), db=5, decode_responses=True)

def connect() -> DBHandler:
    db = DBHandler(auto_commit=True)
    db.connect(
        user=os.environ["POSTGRES_USER"],
        password=os.environ["POSTGRES_PASSWORD"],
        host="postgres",
        port=os.environ["POSTGRES_PORT"],
        db=os.environ["POSTGRES_DB"],
    )
    return db


@celery_app.task(bind=True)
def read_h5ad(self, file_path: str):
    db = connect()
    from . import io
    
    with db as session, StdoutCaptureHandler("general", stdout_redis):
        print("Starting to read h5ad file...")
        io.read_h5ad(db=session, file_path=file_path)
        print("Finished reading h5ad file.")