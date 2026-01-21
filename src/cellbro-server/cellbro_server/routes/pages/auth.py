from fastapi import APIRouter, Cookie, Depends

from cellbro_db import models
from cellbro_db.core.session import AsyncSession

from ...core import exceptions as exc
from ...core.cache import session_cache
from ...core.responses import html_response
from ...core.dependencies import db_session

router = APIRouter(prefix="/auth", tags=["auth", "view"])

@router.get("/login")
async def login(sid: str | None = Cookie(default=None, alias="session_id"), db: AsyncSession = Depends(db_session)):
    if sid:
        if (user_id := await session_cache.get(sid)) is not None:
            if (user := await db.get(models.User.Get(id=user_id))) is not None:
                print(user)
                return await html_response(redirect="/")
    return await html_response("views/login.html")