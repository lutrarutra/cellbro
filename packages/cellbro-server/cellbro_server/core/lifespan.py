import os
from contextlib import asynccontextmanager
from fastapi import FastAPI
import sqlalchemy as sa

from cellbro_worker import tasks
from cellbro_db import types

from .config import settings
from .database import db_handler
from .cache import session_cache, flash_cache, worker_redis


@asynccontextmanager
async def lifespan(app: FastAPI):
    from cellbro_worker.broker import task_broker
    await task_broker.startup()
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
    flash_cache.connect(host=settings.REDIS_HOST, port=settings.REDIS_PORT, db=1)
    worker_redis.connect(host=settings.REDIS_HOST, port=settings.REDIS_PORT, db=2)

    paths = [path for path in  os.listdir(settings.DATA_DIR) if path.endswith(".h5ad") or path.endswith(".h5")]

    if len(paths) == 1:
        if not await worker_redis.is_step_completed(types.ChecklistStep.LOAD):
            print(f"Auto-loading dataset from {paths[0]}", flush=True)
            await tasks.io.read_h5ad.kiq(os.path.join(settings.DATA_DIR, paths[0]))

    yield

    await db_handler.close()
    await db_handler._engine.dispose()
    await session_cache.close()
    await flash_cache.close()
    await worker_redis.close()
    await task_broker.shutdown()