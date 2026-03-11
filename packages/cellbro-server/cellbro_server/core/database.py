from cellbro_db import AsyncDBHandler
from .config import settings

db_handler = AsyncDBHandler()

async def init_db():
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
        import sqlalchemy as sa
        await conn.execute(sa.text("CREATE EXTENSION IF NOT EXISTS pg_trgm;"))
        await conn.run_sync(Base.metadata.create_all)

async def close_db():
    await db_handler.close()
    if db_handler._engine:
        await db_handler._engine.dispose()