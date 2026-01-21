import os
from cellbro_db import DBHandler

from .. import celery_app, queues

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
def create_paper_summary(self, paper_id: int) -> None:
    print(f"Starting create_paper_summary task for paper_id: {paper_id}")
    from . import summary
    db = connect()
    with db as session:
        summary.summarize_pdf(session, paper_id)
    db.close_connection()
    print(f"Completed create_paper_summary task for paper_id: {paper_id}")