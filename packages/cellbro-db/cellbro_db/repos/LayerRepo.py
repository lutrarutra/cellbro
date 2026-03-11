import sqlalchemy as sa
import sqlalchemy.sql as sql

from ..models import Layer

class LayerRepo:
    @classmethod
    def Get(cls, id: int | None = None, name: str | None = None) -> sa.Select[tuple[Layer]]:
        q =  sa.select(Layer)
        if id is not None:
            q = q.where(Layer.id == id)
        elif name is not None:
            q = q.where(Layer.name == name)
        else:
            raise ValueError("Either id or name must be provided")
        return q
    
    @classmethod
    def Create(cls, name: str, dtype: str) -> Layer:
        return Layer(name=name, dtype=dtype)
    
    @classmethod
    def Select(
        cls,
        query: sql.Select[tuple[Layer]] = sa.select(Layer),
        name: str | None = None,
    ) -> sa.Select[tuple[Layer]]:
        if name is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(Layer.name, name.lower()).desc()))
        return query
    