from fastapi import APIRouter, Depends
from cellbro_db import AsyncSession
from cellbro_worker import queues

from ...core.dependencies import db_session, get_user
from ...core import responses
from ...core.context import ctx
from ... import components

router = APIRouter(tags=["home", "view"])

@router.get("/")
async def root(db: AsyncSession = Depends(db_session)):
    task_id = queues.read_h5ad("/app/data/pbmc3k.h5ad")
    select_feature_field = components.inputs.SearchSelectField(
        label="Select Feature",
        search_url=ctx.request.url_for("search_features"),
    )

    select_field = components.inputs.SelectInput(
        label="Select Editor",
        options=[
            ("vscode", "VScode"),
            ("vscode-fork", "VScode fork"),
            ("another-vscode-fork", "Another VScode fork"),
        ],
        required=True
    )

    return await responses.html_response("index.html", select_feature_field=select_feature_field, select_field=select_field)