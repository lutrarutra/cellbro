import sqlalchemy as sa
from typing import Sequence

from cellbro_db import models, AsyncSession, types

async def get_num_observations(db: AsyncSession) -> int:
    return await db.count(sa.select(sa.func.count(models.Observation.id)))

async def get_num_features(db: AsyncSession) -> int:
    return await db.count(sa.select(sa.func.count(models.Feature.id)))

async def get_obs_variables(db: AsyncSession) -> Sequence[models.Variable]:
    return await db.find(models.Variable.Select(layer=types.AnnDataLayerType.OBS, limit=None))

async def get_var_variables(db: AsyncSession) -> Sequence[models.Variable]:
    return await db.find(models.Variable.Select(layer=types.AnnDataLayerType.VAR, limit=None))

async def get_layers(db: AsyncSession) -> Sequence[models.Layer]:
    return await db.find(models.Layer.Select(limit=None))

