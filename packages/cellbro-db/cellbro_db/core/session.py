from typing import TypeVar, Any, Sequence, Iterator, AsyncIterator

import sqlalchemy as sa
from sqlalchemy import sql
from sqlalchemy.ext.asyncio import AsyncSession as SQLAlchemyAsyncSession
from sqlalchemy.orm import Session

from ..models.Base import Base
from .exceptions import ObjectNotFound

T = TypeVar("T", bound=Base)

class AsyncSession(SQLAlchemyAsyncSession):
    async def get_all(
        self,
        statement: sa.Select[tuple[T]],
        order_by: sql.expression.UnaryExpression | None = None,
        limit: int | None = 10,
        options: sa.sql.base.ExecutableOption | None = None,
        offset: int | None = None
    ) -> Sequence[T]:
        """Execute a select statement and return a sequence of objects."""
        if order_by is not None:
            statement = statement.order_by(order_by)
        if offset is not None:
            statement = statement.offset(offset)
        if limit is not None:
            statement = statement.limit(limit)
        if options is not None:
            statement = statement.options(options)

        result = await super().execute(statement)
        return result.scalars().all()
    
    async def page(
        self,
        statement: sa.Select[tuple[T]],
        page: int,
        order_by: sql.expression.UnaryExpression | None = None,
        limit: int = 10,
        options: sa.sql.base.ExecutableOption | None = None
    ) -> tuple[Sequence[T], int]:
        """Execute a select statement and return a page of objects."""
        count = await self.count(statement)
        offset = page * limit
        
        if offset >= count:
            return [], count

        return await self.get_all(statement, limit=limit, offset=offset, order_by=order_by, options=options), count

    async def get_one(self, statement: sa.Select[tuple[T]], options: sa.sql.base.ExecutableOption | None = None) -> T | None:
        """Execute a select statement and return a single object or None."""
        if options is not None:
            statement = statement.options(options)
        result = await super().execute(statement.limit(1))
        return result.scalar_one_or_none()
    
    async def get_or_fail(self, statement: sa.Select[tuple[T]], options: sa.sql.base.ExecutableOption | None = None) -> T:
        """Execute a select statement and return a single object or raise KeyError."""
        if options is not None:
            statement = statement.options(options)
        result = await super().execute(statement)
        if (obj := result.scalar_one_or_none()) is None:
            raise ObjectNotFound()
        return obj

    async def count(self, statement: sa.Select[tuple[T]]) -> int:
        """Execute a count statement and return the integer."""
        result = await super().execute(sa.select(sa.func.count()).select_from(statement.order_by(None).subquery()))
        return result.scalar_one() or 0

    async def save(self, model: T, flush: bool = False) -> T:
        """Add a model to the session and optionally commit."""
        self.add(model)
        if flush:
            await self.flush()
        return model
    
    async def delete(self, model: Any, flush: bool = False) -> None:
        """Delete a model from the session and optionally commit."""
        await super().delete(model)
        if flush:
            await self.flush()

    async def __getitem__(self, statement: sa.Select[tuple[T]]) -> T:
        return await self.get_or_fail(statement)
    
    async def iter(
        self, statement: sa.Select[tuple[T]],
        order_by: sql.expression.UnaryExpression | None = None,
        batch_size: int = 100, options: sa.sql.base.ExecutableOption | None = None
    ) -> AsyncIterator[T]:
        """Execute a select statement and return an async iterator of objects."""
        if order_by is not None:
            statement = statement.order_by(order_by)
        if options is not None:
            statement = statement.options(options)
        
        statement = statement.execution_options(stream_results=True)
        result = await super().stream(statement)
        return result.scalars().yield_per(batch_size)
        

class SyncSession(Session):
    def get_all(
        self,
        statement: sa.Select[tuple[T]],
        limit: int | None = 10,
        order_by: sql.expression.UnaryExpression | None = None,
        options: sa.sql.base.ExecutableOption | None = None,
        offset: int | None = None
    ) -> Sequence[T]:
        """Execute a select statement and return a sequence of objects."""
        if offset is not None:
            statement = statement.offset(offset)
        if order_by is not None:
            statement = statement.order_by(order_by)
        if limit is not None:
            statement = statement.limit(limit)
        if options is not None:
            statement = statement.options(options)

        result = super().execute(statement)
        return result.scalars().all()
    
    def page(
        self,
        statement: sa.Select[tuple[T]], page: int,
        limit: int = 10,
        order_by: sql.expression.UnaryExpression | None = None,
        options: sa.sql.base.ExecutableOption | None = None
    ) -> tuple[Sequence[T], int]:
        """Execute a select statement and return a page of objects."""
        count = self.count(statement)
        offset = page * limit
        
        if offset >= count:
            return [], count
        
        return self.get_all(statement, limit=limit, order_by=order_by, offset=offset, options=options), count
    
    def get_one(self, statement: sa.Select[tuple[T]], options: sa.sql.base.ExecutableOption | None = None) -> T | None:
        """Execute a select statement and return a single object or None."""
        if options is not None:
            statement = statement.options(options)
        result = super().execute(statement.limit(1))
        return result.scalar_one_or_none()
    
    def get_or_fail(self, statement: sa.Select[tuple[T]], options: sa.sql.base.ExecutableOption | None = None) -> T:
        """Execute a select statement and return a single object or raise KeyError."""
        if options is not None:
            statement = statement.options(options)
        result = super().execute(statement.limit(1))
        if (obj := result.scalar_one_or_none()) is None:
            raise ObjectNotFound()
        return obj
    
    def count(self, statement: sa.Select[tuple[T]]) -> int:
        """Execute a count statement and return the integer."""
        result = super().execute(sa.select(sa.func.count()).select_from(statement.order_by(None).subquery()))
        return result.scalar_one() or 0
    
    def save(self, model: T, flush: bool = False) -> T:
        """Add a model to the session and optionally commit."""
        self.add(model)
        if flush:
            self.flush()
        return model
    
    def delete(self, model: Base, flush: bool = False) -> None:
        """Delete a model from the session and optionally commit."""
        super().delete(model)
        if flush:
            self.flush()

    def iter(
        self, statement: sa.Select[tuple[T]],
        order_by: sql.expression.UnaryExpression | None = None,
        options: sa.sql.base.ExecutableOption | None = None,
        batch_size: int = 100,
    ) -> Iterator[T]:
        """Execute a select statement and return an iterator of objects."""
        if order_by is not None:
            statement = statement.order_by(order_by)

        if options is not None:
            statement = statement.options(options)
        
        statement = statement.execution_options(stream_results=True)
        results = super().execute(statement)
        return results.scalars().yield_per(batch_size)


    def __getitem__(self, statement: sa.Select[tuple[T]]) -> T:
        return self.get_or_fail(statement)