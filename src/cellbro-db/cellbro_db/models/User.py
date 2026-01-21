from typing import Self

import sqlalchemy as sa
from sqlalchemy.orm import Mapped, mapped_column
from sqlalchemy.ext.hybrid import hybrid_property

from .Base import Base
from ..types import UserType

class User(Base):
    __tablename__ = "cb_user"

    id: Mapped[int] = mapped_column(sa.Integer, primary_key=True)
    first_name: Mapped[str] = mapped_column(sa.String(64), nullable=False)
    last_name: Mapped[str] = mapped_column(sa.String(64), nullable=False)
    email: Mapped[str] = mapped_column(sa.String(128), nullable=False, unique=True, index=True)
    password: Mapped[str] = mapped_column(sa.String(128), nullable=False)

    type: Mapped[UserType] = mapped_column(sa.Integer, nullable=False, default=UserType.REGULAR)

    @hybrid_property
    def name(self) -> str:  # type: ignore[override]
        return self.first_name + " " + self.last_name
    
    @name.expression
    def name(cls) -> sa.ScalarSelect[str]:
        return sa.select(
            (cls.first_name + " " + cls.last_name)  # type: ignore[arg-type]
        ).correlate(cls).scalar_subquery()
    
    @classmethod
    def Get(cls, id: int | None = None, email: str | None = None) -> sa.Select[tuple[Self]]:
        q =  sa.select(cls)
        if id is not None:
            q = q.where(cls.id == id)
        if email is not None:
            q = q.where(cls.email == email.lower())
        return q
    
    @classmethod
    def Create(cls, first_name: str, last_name: str, email: str, hashed_password: str, type: UserType = UserType.REGULAR) -> "User":
        return cls(
            first_name=first_name,
            last_name=last_name,
            email=email.lower(),
            password=hashed_password,
            type=type,
        )
    
    @classmethod
    def Select(
        cls,
        name: str | None = None,
        limit: int | None = 10, offset: int | None = None,
        sort_by: str | None = None, descending: bool = False,
        page: int | None = None,
    ) -> sa.Select[tuple[Self]]:
        query = sa.select(cls)

        if sort_by is not None:
            attr = getattr(User, sort_by)
            if descending:
                attr = attr.desc()
            query = query.order_by(attr)

        if name is not None:
            query = query.order_by(sa.nulls_last(sa.func.similarity(User.first_name + ' ' + User.last_name, name).desc()))
        

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
            "trgm_cb_user_email_idx",
            sa.text("lower(email) gin_trgm_ops"),
            postgresql_using="gin",
        ),
        sa.Index(
            "trgm_cb_user_first_name_idx",
            sa.text("lower(first_name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
        sa.Index(
            "trgm_cb_user_last_name_idx",
            sa.text("lower(last_name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
        sa.Index(
            "trgm_cb_user_full_name_idx",
            sa.text("lower(first_name || ' ' || last_name) gin_trgm_ops"),
            postgresql_using="gin",
        )
    )