from fastapi import APIRouter, Depends
from cellbro_db import models, AsyncSession

from ...core.dependencies import db_session, get_user
from ...core import responses
from ...core.context import ctx

router = APIRouter(tags=["home", "view"])

@router.get("/")
async def root(db: AsyncSession = Depends(db_session)):
    return await responses.html_response("index.html")