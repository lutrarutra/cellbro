from fastapi import APIRouter, Depends

from ...core import dependencies
from ...core import responses
from ...core.context import ctx
from ... import components
from ...core.cache import worker_redis

router = APIRouter(tags=["qc", "view"])

@router.get("/qc")
async def qc_page():
    return await responses.html_response("base.html", page="qc", page_content_url=ctx.url_for("qc_view"))


@router.get("/views/qc")
async def qc_view(qc_required = Depends(dependencies.qc)):
    return await responses.htmx_response("views/qc.html")
    