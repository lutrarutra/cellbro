import sqlalchemy as sa
import sqlalchemy.sql as sql

from ..types import SoftwareType

from ..models import Software, User

class SoftwareRepo:
    @classmethod
    def Get(cls, id: int | None = None, name: str | None = None) -> sa.Select[tuple[Software]]:
        q =  sa.select(Software)
        if id is not None:
            q = q.where(Software.id == id)
        if name is not None:
            q = q.where(Software.name == name)
        return q
    
    @classmethod
    def Create(cls, name: str, publisher: User, type: SoftwareType) -> Software:
        return Software(
            name=name, publisher=publisher, type=type,
        )
    
    @classmethod
    def Select(
        cls,
        query: sql.Select[tuple[Software]] = sa.select(Software),
        name: str | None = None,
    ) -> sa.Select[tuple[Software]]:
        if name is not None:
            query = query.order_by(sa.func.similarity(Software.name, name).desc())
        return query
