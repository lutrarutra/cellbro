import sqlalchemy as sa
from typing import Sequence

from cellbro_db import models, AsyncSession, types, queries

async def get_num_observations(db: AsyncSession) -> int:
    return await db.count(queries.observation.select())

async def get_num_features(db: AsyncSession) -> int:
    return await db.count(queries.feature.select())

async def get_obs_variables(db: AsyncSession) -> Sequence[models.Variable]:
    return await db.get_all(queries.variable.select(layer=types.AnnDataLayerType.OBS), limit=None, order_by=sa.nulls_last(models.Variable.id.asc()))

async def get_var_variables(db: AsyncSession) -> Sequence[models.Variable]:
    return await db.get_all(queries.variable.select(layer=types.AnnDataLayerType.VAR), limit=None, order_by=models.Variable.name.desc())

async def get_layers(db: AsyncSession) -> Sequence[models.Layer]:
    return await db.get_all(queries.layer.select(), limit=None)

