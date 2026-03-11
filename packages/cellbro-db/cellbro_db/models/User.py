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

    type: Mapped[UserType] = mapped_column(sa.Enum(UserType), nullable=False, default=UserType.REGULAR)

    @hybrid_property
    def name(self) -> str:  # type: ignore[override]
        return self.first_name + " " + self.last_name
    
    @name.expression
    def name(cls) -> sa.ScalarSelect[str]:
        return sa.select(
            (cls.first_name + " " + cls.last_name)  # type: ignore[arg-type]
        ).correlate(cls).scalar_subquery()

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