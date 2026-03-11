import sqlalchemy as sa
import sqlalchemy.sql as sql

from ..models import Observation

class ObservationRepo:
    @classmethod
    def Get(cls, id: int | None = None, name: str | None = None) -> sa.Select[tuple[Observation]]:
        q =  sa.select(Observation)
        if id is not None:
            q = q.where(Observation.id == id)
        elif name is not None:
            q = q.where(Observation.name == name)
        else:
            raise ValueError("Either id or name must be provided")
        return q
    
    @classmethod
    def Create(cls, name: str) -> Observation:
        return Observation(name=name)

    @classmethod
    def Select(
        cls,
        query: sql.Select[tuple[Observation]] = sa.select(Observation),
        name: str | None = None,
    ) -> sa.Select[tuple[Observation]]:
        if name is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(Observation.name, name.lower()).desc()))
        return query