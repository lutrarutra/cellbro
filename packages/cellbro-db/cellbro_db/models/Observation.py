import sqlalchemy as sa
from sqlalchemy.orm import Mapped, mapped_column

from .Base import Base

class Observation(Base):
    __tablename__ = "observation"

    id: Mapped[int] = mapped_column(sa.Integer, primary_key=True)
    name: Mapped[str] = mapped_column(sa.String(128), nullable=False, unique=True, index=True)
    
    __table_args__ = (
        sa.Index(
            "trgm_observation_name_idx",
            sa.text("lower(name) gin_trgm_ops"),
            postgresql_using="gin",
        ),
    )
    
