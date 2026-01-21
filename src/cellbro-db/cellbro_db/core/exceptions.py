class DBException(Exception):
    """Base exception for database errors."""
    def __init__(self, message: str = "Database Error") -> None:
        super().__init__(message)


class ObjectNotFound(DBException):
    """Raised when a requested object is not found in the database."""
    def __init__(self, message: str = "Object not found in the database") -> None:
        super().__init__(message)

