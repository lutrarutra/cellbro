from fastapi import Request, Response, status, HTTPException

from ..core.responses import htmx_response, html_response
from ..core.cache import flash_cache

async def not_authenticated_handler(request: Request, _: Exception) -> Response:
    if request.headers.get("HX-Request"):
        return await htmx_response(redirect="/auth/login")
    return await html_response(redirect="/auth/login")


async def generic_exception_handler(request: Request, exc: Exception) -> Response:
    if request.headers.get("HX-Request"):
        if (sid := request.cookies.get("session_id")):
            await flash_cache.add(sid, category="error", message=f"An internal server error occurred: {str(exc)}")
        return await htmx_response(status=status.HTTP_204_NO_CONTENT)
    return Response(
        content=f"Internal server error: {str(exc)}",
        status_code=exc.status_code if isinstance(exc, HTTPException) else status.HTTP_500_INTERNAL_SERVER_ERROR
    )