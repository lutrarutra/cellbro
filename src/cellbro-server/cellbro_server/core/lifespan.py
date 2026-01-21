from contextlib import asynccontextmanager

from fastapi import FastAPI
import sqlalchemy as sa

from .config import settings
from .database import db_handler
from .cache import session_cache, flash_cache

@asynccontextmanager
async def lifespan(app: FastAPI):
    await db_handler.connect(
        user=settings.POSTGRES_USER,
        password=settings.POSTGRES_PASSWORD,
        host=settings.POSTGRES_HOST,
        db=settings.POSTGRES_DB,
        port=settings.POSTGRES_PORT
    )
    if db_handler._engine is None:
        raise Exception("DB connection could not be established")
    
    async with db_handler._engine.begin() as conn:
        from cellbro_db.models.Base import Base
        await conn.execute(sa.text("CREATE EXTENSION IF NOT EXISTS pg_trgm;"))
        await conn.run_sync(Base.metadata.create_all)

    session_cache.connect(host=settings.REDIS_HOST, port=settings.REDIS_PORT, db=0)
    await flash_cache.connect(host=settings.REDIS_HOST, port=settings.REDIS_PORT, db=1)
    
    yield
    await db_handler.close()
    await db_handler._engine.dispose()
    await session_cache.close()
    await flash_cache.close()