from fastapi import Depends, APIRouter

router = APIRouter(prefix="/htmx/search", tags=["search", "htmx"])

from cellbro_db import models
from cellbro_db.core.session import AsyncSession

from ...core.context import ctx
from ...core import responses
from ... import forms, logic
from ...core.dependencies import db_session, get_user


@router.get("/search_features")
async def search_features(word: str, db: AsyncSession = Depends(db_session)):
    results = await db.find(models.Feature.Select(word=word, sort_by=models.Feature.name if not word else None))
    return await responses.htmx_response("components/search/results.html", results=results)