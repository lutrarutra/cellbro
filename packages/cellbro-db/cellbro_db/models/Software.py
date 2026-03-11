from typing import TYPE_CHECKING

import sqlalchemy as sa
import sqlalchemy.orm as orm
from sqlalchemy.orm import Mapped, mapped_column

from .Base import Base
from ..types import SoftwareType


if TYPE_CHECKING:
    from .User import User


class Software(Base):
    __tablename__ = "software"

    id: Mapped[int] = mapped_column(sa.Integer, primary_key=True)
    name: Mapped[str] = mapped_column(sa.String(128), nullable=False, unique=True, index=True)

    publisher_id: Mapped[int] = mapped_column(sa.Integer, sa.ForeignKey("cb_user.id"), nullable=False)
    publisher: Mapped["User"] = orm.relationship("User", backref="published_software")

    type: Mapped[SoftwareType] = mapped_column(sa.Enum(SoftwareType), nullable=False)

    __table_args__ = (
        sa.Index(
            "trgm_software_name_idx",
            sa.text("lower(name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
    )