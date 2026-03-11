from typing import Self

import sqlalchemy as sa
import sqlalchemy.orm as orm
from sqlalchemy.orm import Mapped, mapped_column

from .Base import Base

class Layer(Base):
    __tablename__ = "layer"

    id: Mapped[int] = mapped_column(sa.Integer, primary_key=True)
    name: Mapped[str] = mapped_column(sa.String(128), nullable=False, unique=True, index=True)
    dtype: Mapped[str] = mapped_column(sa.String(64), nullable=True, unique=False, index=False)
    
    __table_args__ = (
        sa.Index(
            "trgm_layer_name_idx",
            sa.text("lower(name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
    )