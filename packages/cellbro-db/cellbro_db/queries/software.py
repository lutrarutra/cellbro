import sqlalchemy as sa
import sqlalchemy.sql as sql

from ..types import SoftwareType

from ..models import Software, User

def get(id: int | None = None, name: str | None = None) -> sa.Select[tuple[Software]]:
    q =  sa.select(Software)
    if id is not None:
        q = q.where(Software.id == id)
    if name is not None:
        q = q.where(Software.name == name)
    return q

def create(name: str, publisher: User, type: SoftwareType) -> Software:
    return Software(
        name=name, publisher=publisher, type=type,
    )

def select(
    name: str | None = None,
    type: SoftwareType | None = None,
    type_in: list[SoftwareType] | None = None,
    query: sql.Select[tuple[Software]] = sa.select(Software),
) -> sa.Select[tuple[Software]]:
    if type is not None:
        query = query.where(Software.type == type)
    if type_in is not None:
        query = query.where(Software.type.in_(type_in))
    if name is not None:
        query = query.order_by(sa.func.similarity(Software.name, name).desc())
    return query
