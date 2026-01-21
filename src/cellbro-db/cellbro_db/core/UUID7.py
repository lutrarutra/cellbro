from sqlalchemy.dialects.postgresql import UUID as PG_UUID
from sqlalchemy.types import TypeDecorator

class UUID7(TypeDecorator):
    """
    Postgres-compatible UUID7 type.
    Generates a time-ordered UUID7 by default.
    """
    impl = PG_UUID
    cache_ok = True

    def process_bind_param(self, value, dialect):
        if value is None:
            return None
        if isinstance(value, str):
            return value
        return str(value)