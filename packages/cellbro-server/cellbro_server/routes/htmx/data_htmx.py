from fastapi import Depends, APIRouter

from cellbro_db.core.session import AsyncSession
from cellbro_db import queries

from ...core import responses
from ... import logic
from ...core.dependencies import db_session, dataset

router = APIRouter(prefix="/htmx/data", tags=["data", "htmx"])

@router.get("/obs/summary")
async def get_obs_summary(db: AsyncSession = Depends(db_session), _ = Depends(dataset)):
    users = await db.get_all(queries.user.select())
    variables = await logic.dataset.get_obs_variables(db)
    print(f"variables: {variables}")
    return await responses.htmx_response("components/mini/variable-summary.html", variables=variables)


@router.get("/var/summary")
async def get_var_summary(db: AsyncSession = Depends(db_session), _ = Depends(dataset)):
    variables = await logic.dataset.get_var_variables(db)
    print(f"variables: {variables}")
    return await responses.htmx_response("components/mini/variable-summary.html", variables=variables)


@router.get("/layers/summary")
async def get_layer_summary(db: AsyncSession = Depends(db_session), _ = Depends(dataset)):
    layers = await logic.dataset.get_layers(db)
    return await responses.htmx_response("components/mini/layer-summary.html", layers=layers)

@router.get("/dataset/summary")
async def get_dataset_summary(db: AsyncSession = Depends(db_session), _ = Depends(dataset)):
    num_observations = await logic.dataset.get_num_observations(db)
    num_features = await logic.dataset.get_num_features(db)
    print(f"num_observations: {num_observations}, num_features: {num_features}")
    count = 0
    async for feature in await db.iter(queries.feature.select()):
        count += 1
    print(f"count: {count}")
    return await responses.htmx_response("components/mini/dataset-summary.html", num_observations=num_observations, num_features=num_features)