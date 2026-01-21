import os
from cellbro_db import DBHandler

from .. import celery_app

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
    from . import io
    db = connect()
    with db as session:
        io.read_h5ad(db=session, file_path=file_path)