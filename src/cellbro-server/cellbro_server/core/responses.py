from typing import Any

from fastapi.responses import HTMLResponse, RedirectResponse, Response

from ..core.templates import render_template
from ..core.context import ctx

async def get_request_context() -> dict[str, Any]:
    request = ctx.request
    return {
        "request": request,
        "current_user": await ctx.current_user,
    }

async def html_response(
    template: str | None = None, 
    redirect: str | None = None, 
    status: int = 200, 
    **context
) -> Response:
    if redirect:
        return RedirectResponse(url=redirect, status_code=303)
    
    content = ""
    if template is not None:
        content = await render_template(template, **context | await get_request_context())
    
    return HTMLResponse(
        content=content,
        status_code=status,
        headers={"Content-Type": "text/html; charset=utf-8"}
    )

async def htmx_response(
    template: str | None = None, 
    status: int = 200, 
    redirect: str | None = None, 
    re_target: str | None = None, 
    re_swap: str | None = None, # Added for completeness
    **context
) -> Response:
    headers = {"HX-Trigger": "contentUpdated"}
    
    if redirect:
        headers["HX-Redirect"] = redirect
        return HTMLResponse(status_code=204, headers=headers)

    if re_target:
        headers["HX-Retarget"] = re_target
    
    if re_swap:
        headers["HX-Reswap"] = re_swap
    
    content = ""
    if template is not None:
        content = await render_template(template, **context | await get_request_context())
    elif status == 200:
        status = 204

    return HTMLResponse(
        content=content,
        status_code=status,
        headers=headers
    )