from typing import Sequence, Callable

import sqlalchemy as sa
from sqlalchemy.orm import Session
from sqlalchemy.ext.asyncio import AsyncSession
from sqlalchemy.orm import Query

from ...types import UserType
from ... import models

class UserBP:
    def __init__(self, session: Session | AsyncSession):
        self.session = session

    @classmethod
    def where(
        cls,
        query: Query,
        type: UserType | None = None,
        type_in: Sequence[UserType] | None = None,
        custom_query: Callable[[Query], Query] | None = None,
    ):
        if type is not None:
            query = query.filter(models.User.type == type)
        if type_in is not None:
            query = query.filter(models.User.type.in_(type_in))
        if custom_query is not None:
            query = custom_query(query)
        return query
    
    def create_user(
        self,
        first_name: str,
        last_name: str,
        email: str,
        hashed_password: str,
        type: UserType = UserType.REGULAR,
    ) -> models.User:
        new_user = models.User.create(
            first_name=first_name,
            last_name=last_name,
            email=email,
            hashed_password=hashed_password,
            type=type,
        )
        self.session.add(new_user)
        return new_user
        