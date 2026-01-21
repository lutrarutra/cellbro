from uuid import uuid4
from fastapi import Request, Cookie, Depends

from cellbro_db.core.session import AsyncSession
from cellbro_db import models

from .database import db_handler
from .cache import session_cache
from . import exceptions as exc

async def db_session(request: Request):
    async with db_handler.get_session() as session:
        request.state.db = session
        yield session
        if bool(session.dirty) or bool(session.new) or bool(session.deleted):
            try:
                await session.commit()
            except Exception as e:
                await session.rollback()
                raise e

async def get_sid(sid: str | None = Cookie(default=None, alias="session_id")) -> str:
    if not sid:
        sid = uuid4().hex
    return sid

async def get_user(db: AsyncSession = Depends(db_session), sid: str | None = Cookie(default=None, alias="session_id")) -> models.User:
    if not sid:
        raise exc.NotAuthenticatedException()

    user_id = await session_cache.get(sid)
    if not user_id:
        raise exc.NotAuthenticatedException()

    if not (user := await db.get(models.User.Get(id=user_id))):
        raise exc.NotAuthenticatedException()
    return user