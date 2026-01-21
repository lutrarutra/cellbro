from typing import TypeVar, Any, Sequence
import sqlalchemy as sa
from sqlalchemy.ext.asyncio import AsyncSession as SQLAlchemyAsyncSession
from sqlalchemy.orm import Session

from ..models.Base import Base
from .exceptions import ObjectNotFound

T = TypeVar("T", bound=Base)

class AsyncSession(SQLAlchemyAsyncSession):
    async def find(self, statement: sa.Select[tuple[T]]) -> Sequence[T]:
        """Execute a select statement and return a sequence of objects."""
        result = await super().execute(statement)
        return result.scalars().all()

    async def get(self, statement: sa.Select[tuple[T]]) -> T | None:
        """Execute a select statement and return a single object or None."""
        result = await super().execute(statement)
        return result.scalar_one_or_none()
    
    async def get_or_fail(self, statement: sa.Select[tuple[T]]) -> T:
        """Execute a select statement and return a single object or raise KeyError."""
        result = await super().execute(statement)
        obj = result.scalar_one_or_none()
        if obj is None:
            raise ObjectNotFound()
        return obj

    async def count(self, statement: sa.Select[tuple[int]]) -> int:
        """Execute a count statement and return the integer."""
        result = await super().execute(statement)
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

class SyncSession(Session):
    def find(self, statement: sa.Select[tuple[T]]) -> Sequence[T]:
        """Execute a select statement and return a sequence of objects."""
        result = super().execute(statement)
        return result.scalars().all()
    
    def get(self, statement: sa.Select[tuple[T]]) -> T | None:
        """Execute a select statement and return a single object or None."""
        result = super().execute(statement)
        return result.scalar_one_or_none()
    
    def get_or_fail(self, statement: sa.Select[tuple[T]]) -> T:
        """Execute a select statement and return a single object or raise KeyError."""
        result = super().execute(statement)
        obj = result.scalar_one_or_none()
        if obj is None:
            raise ObjectNotFound()
        return obj
    
    def count(self, statement: sa.Select[tuple[int]]) -> int:
        """Execute a count statement and return the integer."""
        result = super().execute(statement)
        return result.scalar_one() or 0
    
    def save(self, model: T, flush: bool = False) -> T:
        """Add a model to the session and optionally commit."""
        self.add(model)
        if flush:
            self.flush()
        return model
    
    def delete(self, model: Base, flush: bool = False) -> None:
        """Delete a model from the session and optionally commit."""
        self.delete(model)
        if flush:
            self.flush()

    def __getitem__(self, statement: sa.Select[tuple[T]]) -> T:
        return self.get_or_fail(statement)