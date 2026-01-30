from fastapi import APIRouter

from ...core import responses
from ...core.context import ctx

router = APIRouter(tags=["dashboard", "view"])

@router.get("/")
async def dashboard():
    return await responses.html_response("base.html", page="dashboard", page_content_url=ctx.url_for("dashboard_view"))

@router.get("/views/dashboard")
async def dashboard_view():
    return await responses.htmx_response("views/dashboard.html")