import sqlalchemy as sa
from sqlalchemy.orm import DeclarativeBase

class Base(DeclarativeBase):    
    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(id={getattr(self, 'id', None)})"

    def __str__(self) -> str:
        return self.__repr__()