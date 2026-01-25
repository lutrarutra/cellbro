from fastapi import APIRouter, Depends

from ...core import dependencies
from ...core import responses
from ...core.context import ctx
from ... import components
from ...core.cache import worker_redis

router = APIRouter(tags=["qc", "view"])

@router.get("/qc")
async def qc_page(qc_required = Depends(dependencies.qc)):
    return await responses.html_response("views/qc.html")
    