from fastapi import Depends, APIRouter
import sqlalchemy as sa

router = APIRouter(prefix="/htmx/data", tags=["data", "htmx"])

from cellbro_db import models
from cellbro_db.core.session import AsyncSession

from ...core.context import ctx
from ...core import responses
from ... import forms, logic
from ...core.dependencies import db_session


@router.get("/get")
async def get_dataset_status(db: AsyncSession = Depends(db_session)):
    num_features = await db.execute(sa.select(sa.func.count(models.Feature.id)))
    num_observations = await db.execute(sa.select(sa.func.count(models.Observation.id)))
    return await responses.htmx_response("components/mini/dataset-status.html", num_features=num_features.scalar(), num_observations=num_observations.scalar())