from fastapi import Depends, APIRouter, Request

router = APIRouter(prefix="/htmx/users", tags=["users", "htmx"])

from ...core.templates import j2


@router.get("/create")
async def create(request: Request):
    return j2.TemplateResponse("forms/models/user.html", {"request": request})

@router.post("/create")
async def create_post(request: Request):
    pass
    
