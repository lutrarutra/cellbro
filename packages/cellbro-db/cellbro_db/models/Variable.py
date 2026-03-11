import sqlalchemy as sa
from sqlalchemy.orm import Mapped, mapped_column

from .Base import Base
from ..types import AnnDataLayerType, VariableType

class Variable(Base):
    __tablename__ = "variable"

    id: Mapped[int] = mapped_column(sa.Integer, primary_key=True)
    name: Mapped[str] = mapped_column(sa.String(128), nullable=False, unique=False, index=True)
    
    layer: Mapped[AnnDataLayerType] = mapped_column(sa.String(16), nullable=False, unique=False, index=True)
    type: Mapped[VariableType] = mapped_column(sa.Enum(VariableType, native_enum=False), nullable=False, index=True)
    
    __table_args__ = (
        sa.Index(
            "trgm_variable_name_idx",
            sa.text("lower(name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
        sa.UniqueConstraint("name", "layer", name="uq_variable_name_layer"),
    )