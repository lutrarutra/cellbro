import sqlalchemy as sa
from typing import Sequence

from cellbro_db import models, AsyncSession, types, repos

async def get_num_observations(db: AsyncSession) -> int:
    return await db.count(repos.ObservationRepo.Select())

async def get_num_features(db: AsyncSession) -> int:
    return await db.count(repos.FeatureRepo.Select())

async def get_obs_variables(db: AsyncSession) -> Sequence[models.Variable]:
    return await db.get_all(repos.VariableRepo.Select(layer=types.AnnDataLayerType.OBS), limit=None, order_by=sa.nulls_last(models.Variable.id.asc()))

async def get_var_variables(db: AsyncSession) -> Sequence[models.Variable]:
    return await db.get_all(repos.VariableRepo.Select(layer=types.AnnDataLayerType.VAR), limit=None, order_by=models.Variable.name.desc())

async def get_layers(db: AsyncSession) -> Sequence[models.Layer]:
    return await db.get_all(repos.LayerRepo.Select(), limit=None)

