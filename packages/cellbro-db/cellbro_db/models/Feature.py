import sqlalchemy as sa
from sqlalchemy.orm import Mapped, mapped_column

from .Base import Base

class Feature(Base):
    __tablename__ = "feature"

    id: Mapped[int] = mapped_column(sa.Integer, primary_key=True)
    identifier: Mapped[str] = mapped_column(sa.String(128), nullable=False, unique=True, index=True)
    name: Mapped[str] = mapped_column(sa.String(128), nullable=False)
    
    __table_args__ = (
        sa.Index(
            "trgm_feature_name_idx",
            sa.text("lower(name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
        sa.Index(
            "trgm_feature_identifier_idx",
            sa.text("lower(identifier) gin_trgm_ops"),
            postgresql_using="gin",
        ),
        sa.Index(
            "trgm_feature_name_identifier_idx",
            sa.text("lower(identifier || ' ' || name) gin_trgm_ops"),
            postgresql_using="gin",
        )
    )
    
