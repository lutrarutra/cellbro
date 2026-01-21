from fastapi import HTTPException, status

class NotAuthenticatedException(HTTPException):
    def __init__(self, detail: str = "Not authenticated", headers: dict | None = None):
        super().__init__(
            status_code=status.HTTP_401_UNAUTHORIZED, 
            detail=detail,
            headers=headers
        )

class ItemNotFoundException(HTTPException):
    def __init__(self, item_id: str):
        super().__init__(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Item with ID {item_id} was not found"
        )

class ResourceNotFoundException(HTTPException):
    def __init__(self, resource_name: str):
        super().__init__(
            status_code=status.HTTP_404_NOT_FOUND,
            detail=f"Resource '{resource_name}' was not found"
        )

class PermissionDeniedException(HTTPException):
    def __init__(self, detail: str = "Permission denied"):
        super().__init__(
            status_code=status.HTTP_403_FORBIDDEN,
            detail=detail
        )

class InvalidCredentialsException(Exception):
    pass