from uuid import uuid4
from fastapi import Request, Cookie, Depends

from cellbro_db import models, types, AsyncSession, queries

from .database import db_handler
from .cache import session_cache, worker_redis
from . import exceptions as exc

async def dataset(request: Request):
    if not await worker_redis.is_step_completed(types.ChecklistStep.LOAD):
        pass
    
async def qc(request: Request):
    if not await worker_redis.is_step_completed(types.ChecklistStep.LOAD):
        raise exc.DatasetNotLoadedException()
    if not await worker_redis.is_step_completed(types.ChecklistStep.QC):
        raise exc.QCNotCompletedException()

async def db_session(request: Request):
    async with db_handler.get_session() as session:
        request.state.db = session
        try:
            yield session
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

    if not (user := await db.get_one(queries.user.get(id=user_id))):
        raise exc.NotAuthenticatedException()
    return user