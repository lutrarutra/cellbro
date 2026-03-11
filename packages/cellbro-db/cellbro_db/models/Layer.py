from typing import Self

import sqlalchemy as sa
import sqlalchemy.orm as orm
from sqlalchemy.orm import Mapped, mapped_column

from .Base import Base

class Layer(Base):
    __tablename__ = "layer"

    id: Mapped[int] = mapped_column(sa.Integer, primary_key=True)
    name: Mapped[str] = mapped_column(sa.String(128), nullable=False, unique=True, index=True)
    dtype: Mapped[str] = mapped_column(sa.String(64), nullable=True, unique=False, index=False)

    @classmethod
    def Get(cls, id: int | None = None, name: str | None = None) -> sa.Select[tuple[Self]]:
        q =  sa.select(cls)
        if id is not None:
            q = q.where(cls.id == id)
        elif name is not None:
            q = q.where(cls.name == name)
        else:
            raise ValueError("Either id or name must be provided")
        return q
    
    @classmethod
    def Create(cls, name: str, dtype: str) -> "Layer":
        return cls(name=name, dtype=dtype)
    
    @classmethod
    def Select(
        cls,
        name: str | None = None,
        limit: int | None = 10, offset: int | None = None,
        sort_by: orm.InstrumentedAttribute | str | None = None, descending: bool = False,
        page: int | None = None,
    ) -> sa.Select[tuple[Self]]:
        query = sa.select(cls)

        if sort_by is not None:
            if isinstance(sort_by, str):
                attr = getattr(cls, sort_by)
            else:
                attr = sort_by
            if descending:
                attr = attr.desc()
            query = query.order_by(attr)

        if name is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(cls.name, name.lower()).desc()))
        
        if page is not None:
            if limit is None:
                raise ValueError("Limit must be set if page is set.")
            offset = page * limit
        
        if offset is not None:
            query = query.offset(offset)
        
        if limit is not None:
            query = query.limit(limit)
        
        return query
    
    @classmethod
    def Page(
        cls,
        name: str | None = None,
        limit: int | None = 10,
        sort_by: str | None = None,
        descending: bool = False,
        page: int | None = None,
    ) -> tuple[sa.Select[tuple[Self]], sa.Select[tuple[int]]]:
        
        query = cls.select(
            name=name,
            sort_by=sort_by,
            descending=descending,
            limit=None,
            offset=None,
            page=None
        )

        counts_query = sa.select(sa.func.count()).select_from(query.subquery())

        if page is not None and limit is not None:
            query = query.offset(page * limit)
        if limit is not None:
            query = query.limit(limit)

        return query, counts_query
    
    __table_args__ = (
        sa.Index(
            "trgm_layer_name_idx",
            sa.text("lower(name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
    )