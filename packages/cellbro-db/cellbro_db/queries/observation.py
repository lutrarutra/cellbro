import sqlalchemy as sa
import sqlalchemy.sql as sql

from ..models import Observation

def get(id: int | None = None, name: str | None = None) -> sa.Select[tuple[Observation]]:
    q =  sa.select(Observation)
    if id is not None:
        q = q.where(Observation.id == id)
    elif name is not None:
        q = q.where(Observation.name == name)
    else:
        raise ValueError("Either id or name must be provided")
    return q

def create(name: str) -> Observation:
    return Observation(name=name)

def select(
    query: sql.Select[tuple[Observation]] = sa.select(Observation),
    name: str | None = None,
) -> sa.Select[tuple[Observation]]:
    if name is not None:
        query = query.order_by(sa.nulls_last(sa.func.similarity(Observation.name, name.lower()).desc()))
    return query