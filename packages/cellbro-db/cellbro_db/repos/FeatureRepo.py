import sqlalchemy as sa
import sqlalchemy.sql as sql

from ..models import Feature

class FeatureRepo:
    @classmethod
    def Get(cls, id: int | None = None, identifier: str | None = None) -> sa.Select[tuple[Feature]]:
        q =  sa.select(Feature)
        if id is not None:
            q = q.where(Feature.id == id)
        elif identifier is not None:
            q = q.where(Feature.identifier == identifier)
        else:
            raise ValueError("Either id or identifier must be provided")
        return q
    
    @classmethod
    def Create(cls, identifier: str, name: str) -> Feature:
        return Feature(identifier=identifier, name=name)

    @classmethod
    def Select(
        cls,
        word: str | None = None,
        name: str | None = None,
        identifier: str | None = None,
        sort_by: sql.expression.UnaryExpression | None = None,
    ) -> sa.Select[tuple[Feature]]:
        query = sa.select(Feature)

        if name is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(Feature.name, name.lower()).desc()))
        elif identifier is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(Feature.identifier, identifier.lower()).desc()))
        elif word is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(Feature.identifier + ' ' + Feature.name, word.lower()).desc()))
        elif sort_by is not None:
            query = query.order_by(sort_by)

        return query