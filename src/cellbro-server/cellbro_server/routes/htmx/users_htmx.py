from fastapi import Depends, APIRouter, Request

router = APIRouter(prefix="/htmx/users", tags=["users", "htmx"])

from ...core.responses import htmx_response


@router.get("/create")
async def create(request: Request):
    return await htmx_response("forms/models/user.html", request=request)

@router.post("/create")
async def create_post(request: Request):
    pass
    
