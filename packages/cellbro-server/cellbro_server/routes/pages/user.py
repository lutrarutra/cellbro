from fastapi import Depends, APIRouter

from cellbro_db.core.session import AsyncSession
from cellbro_db import models

from ...core import exceptions as exc
from ...core.dependencies import db_session
from ...core.responses import html_response

router = APIRouter(prefix="/users", tags=["users", "view"])

@router.get("/")
async def users(db: AsyncSession = Depends(db_session)):
    users = await db.find(models.User.Select())
    return await html_response("views/users.html", users=users)

@router.get("/{user_id}")
async def user(user_id: int, db: AsyncSession = Depends(db_session)):
    if (user := await db.get_one(models.User.Get(user_id))) is None:
        raise exc.ItemNotFoundException(f"User with ID '{user_id}' not found.")
    return await html_response("views/user.html", user=user)