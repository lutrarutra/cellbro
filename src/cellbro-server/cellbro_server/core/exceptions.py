from fastapi import HTTPException, status, Request, Response

from .context import ctx
from .responses import htmx_response, html_response
from .cache import flash_cache
from .. import forms
from ..core import responses

async def generic_exception_handler(request: Request, exc: Exception) -> Response:
    if request.headers.get("HX-Request"):
        if (sid := request.cookies.get("session_id")):
            await flash_cache.add(sid, category="error", message=f"An internal server error occurred: {str(exc)}")
        return await htmx_response(status=status.HTTP_204_NO_CONTENT)
    return Response(
        content=f"Internal server error: {str(exc)}",
        status_code=exc.status_code if isinstance(exc, HTTPException) else status.HTTP_500_INTERNAL_SERVER_ERROR
    )


class CellBroServerException(Exception):
    @staticmethod
    async def handler(request: Request, _: Exception) -> Response:
        return await generic_exception_handler(request, _)



class NotAuthenticatedException(HTTPException):
    def __init__(self, detail: str = "Not authenticated", headers: dict | None = None):
        super().__init__(
            status_code=status.HTTP_401_UNAUTHORIZED, 
            detail=detail,
            headers=headers
        )

    @staticmethod
    async def handler(request: Request, _: Exception) -> Response:
        if request.headers.get("HX-Request"):
            return await htmx_response(redirect="/auth/login")
        return await html_response(redirect="/auth/login")

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

class InvalidCredentialsException(CellBroServerException):
    pass

class StepRequirementNotMetException(CellBroServerException):
    pass

class DatasetNotLoadedException(StepRequirementNotMetException):
    @staticmethod
    async def handler(request: Request, _: Exception) -> Response:
        if request.headers.get("HX-Request"):
            return await htmx_response(redirect=ctx.request.url_for("dashboard"))
        return await html_response(redirect=ctx.request.url_for("dashboard"))

class QCNotCompletedException(StepRequirementNotMetException):
    @staticmethod
    async def handler(request: Request, _: Exception) -> Response:
        form = forms.steps.QCForm()
        if request.headers.get("HX-Request"):
            return await form.make_response()
        return await responses.html_response("views/qc.html", modal_form=form)

