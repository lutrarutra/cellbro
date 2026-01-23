from typing import Self

import sqlalchemy as sa
import sqlalchemy.orm as orm
from sqlalchemy.orm import Mapped, mapped_column

from .Base import Base
from ..types import AnnDataLayerType, VariableType

class Variable(Base):
    __tablename__ = "variable"

    id: Mapped[int] = mapped_column(sa.Integer, primary_key=True)
    name: Mapped[str] = mapped_column(sa.String(128), nullable=False, unique=False, index=True)
    
    layer: Mapped[AnnDataLayerType] = mapped_column(sa.String(16), nullable=False, unique=False, index=True)
    type: Mapped[VariableType] = mapped_column(sa.Enum(VariableType, native_enum=False), nullable=False, index=True)

    @classmethod
    def Get(cls, id: int) -> sa.Select[tuple[Self]]:
        q =  sa.select(cls)
        q = q.where(cls.id == id)
        return q
    
    @classmethod
    def GetByNameAndLayer(cls, name: str, layer: AnnDataLayerType) -> sa.Select[tuple[Self]]:
        q =  sa.select(cls)
        q = q.where(cls.name == name)
        q = q.where(cls.layer == layer)
        return q
    
    @classmethod
    def Create(cls, name: str, layer: AnnDataLayerType, type: VariableType) -> "Variable":
        return cls(name=name, layer=layer, type=type)
    
    @classmethod
    def Select(
        cls,
        name: str | None = None,
        layer: AnnDataLayerType | None = None,
        type: VariableType | None = None,
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

        if layer is not None:
            query = query.where(cls.layer == layer)

        if type is not None:
            query = query.where(cls.type == type)

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
            "trgm_variable_name_idx",
            sa.text("lower(name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
        sa.UniqueConstraint("name", "layer", name="uq_variable_name_layer"),
    )