import uuid

from cellbro_db import models, AsyncSession, types
from cellbro_db.core.session import AsyncSession

from ..core.cache import session_cache
from ..core import exceptions as exc
from ..core.context import ctx
from ..core.config import settings
from ..core import secrets

async def login(db: AsyncSession, email: str, password: str) -> models.User:
    if (user := await db.get_one(models.User.Get(email=email))) is None:
        raise exc.InvalidCredentialsException("Invalid email or password")
    
    if not secrets.verify_password(password, user.password):
        raise exc.InvalidCredentialsException("Invalid email or password")
    
    session_id = ctx.sid or uuid.uuid4().hex
        
    await session_cache.set(session_id, user.id)
    ctx.response.set_cookie(key="session_id", value=session_id, max_age=settings.SESSION_EXPIRE_SECONDS, httponly=True)

    return user

async def logout(user_id: int) -> None:
    num_sessions_closed = await session_cache.delete_user_sessions(user_id=user_id)
    new_session_id = uuid.uuid4().hex
    ctx.request.cookies["session_id"] = new_session_id
    ctx.response.set_cookie(key="session_id", value=new_session_id, max_age=settings.SESSION_EXPIRE_SECONDS, httponly=True)


async def register(
    db: AsyncSession, email: str, password: str,
    first_name: str, last_name: str,
    account_type: types.UserType
) -> models.User:    
    hashed_password = secrets.hash_password(password)
    user = await db.save(models.User.Create(
        email=email, hashed_password=hashed_password,
        first_name=first_name, last_name=last_name,
        type=account_type
    ))
    return user
    
    