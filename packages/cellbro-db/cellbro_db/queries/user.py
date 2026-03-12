import sqlalchemy as sa
import sqlalchemy.sql as sql

from ..types import UserType
from ..models import User


def get(id: int | None = None, email: str | None = None) -> sa.Select[tuple[User]]:
    q =  sa.select(User)
    if id is not None:
        q = q.where(User.id == id)
    if email is not None:
        q = q.where(User.email == email.lower())
    return q

def create(first_name: str, last_name: str, email: str, hashed_password: str, type: UserType = UserType.REGULAR) -> User:
    return User(
        first_name=first_name,
        last_name=last_name,
        email=email.lower(),
        password=hashed_password,
        type=type,
    )

def select(
    type: UserType | None = None,
    type_in: list[UserType] | None = None,
    query: sql.Select[tuple[User]] = sa.select(User),
    name: str | None = None,
) -> sa.Select[tuple[User]]:
    
    if type is not None:
        query = query.where(User.type == type)
    if type_in is not None:
        query = query.where(User.type.in_(type_in))
    
    if name is not None:
        query = query.order_by(sa.nulls_last(sa.func.similarity(User.first_name + ' ' + User.last_name, name).desc()))
    return query