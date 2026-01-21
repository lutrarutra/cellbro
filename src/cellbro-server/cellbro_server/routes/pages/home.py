from fastapi import APIRouter, Depends
from cellbro_db import AsyncSession
from cellbro_worker import queues

from ...core.dependencies import db_session, get_user
from ...core import responses
from ...core.context import ctx

router = APIRouter(tags=["home", "view"])

@router.get("/")
async def root(db: AsyncSession = Depends(db_session)):
    # task = queues.read_h5ad("/app/data/pbmc3k_raw.h5ad")
    # print(task)
    return await responses.html_response("index.html")