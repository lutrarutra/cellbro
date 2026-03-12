import sqlalchemy as sa
import sqlalchemy.sql as sql

from ..models import Feature

def get(id: int | None = None, identifier: str | None = None) -> sa.Select[tuple[Feature]]:
    q =  sa.select(Feature)
    if id is not None:
        q = q.where(Feature.id == id)
    elif identifier is not None:
        q = q.where(Feature.identifier == identifier)
    else:
        raise ValueError("Either id or identifier must be provided")
    return q

def create(identifier: str, name: str) -> Feature:
    return Feature(identifier=identifier, name=name)

def select(
    word: str | None = None,
    name: str | None = None,
    identifier: str | None = None,
    query: sql.Select[tuple[Feature]] = sa.select(Feature),
) -> sa.Select[tuple[Feature]]:
    if name is not None:
        query = query.order_by(sa.nulls_last(sa.func.similarity(Feature.name, name.lower()).desc()))
    elif identifier is not None:
        query = query.order_by(sa.nulls_last(sa.func.similarity(Feature.identifier, identifier.lower()).desc()))
    elif word is not None:
        query = query.order_by(sa.nulls_last(sa.func.similarity(Feature.identifier + ' ' + Feature.name, word.lower()).desc()))
    return query