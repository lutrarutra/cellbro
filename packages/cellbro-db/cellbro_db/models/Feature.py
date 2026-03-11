from typing import Self

import sqlalchemy as sa
import sqlalchemy.orm as orm
from sqlalchemy.orm import Mapped, mapped_column

from .Base import Base

class Feature(Base):
    __tablename__ = "feature"

    id: Mapped[int] = mapped_column(sa.Integer, primary_key=True)
    identifier: Mapped[str] = mapped_column(sa.String(128), nullable=False, unique=True, index=True)
    name: Mapped[str] = mapped_column(sa.String(128), nullable=False)


    @classmethod
    def Get(cls, id: int | None = None, identifier: str | None = None) -> sa.Select[tuple[Self]]:
        q =  sa.select(cls)
        if id is not None:
            q = q.where(cls.id == id)
        elif identifier is not None:
            q = q.where(cls.identifier == identifier)
        else:
            raise ValueError("Either id or identifier must be provided")
        return q
    
    @classmethod
    def Create(cls, identifier: str, name: str) -> "Feature":
        return cls(identifier=identifier, name=name)
    

    @classmethod
    def Select(
        cls,
        word: str | None = None,
        name: str | None = None,
        identifier: str | None = None,
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
        if identifier is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(cls.identifier, identifier.lower()).desc()))
        if word is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(cls.identifier + ' ' + cls.name, word.lower()).desc()))

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
        sort_by: orm.InstrumentedAttribute | str | None = None,
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
            "trgm_feature_name_idx",
            sa.text("lower(name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
        sa.Index(
            "trgm_feature_identifier_idx",
            sa.text("lower(identifier) gin_trgm_ops"),
            postgresql_using="gin",
        ),
        sa.Index(
            "trgm_feature_name_identifier_idx",
            sa.text("lower(identifier || ' ' || name) gin_trgm_ops"),
            postgresql_using="gin",
        )
    )
    
