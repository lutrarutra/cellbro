from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from .core import middleware as mw, handlers, context
from .core.lifespan import lifespan
from .core.exceptions import NotAuthenticatedException
from . import routes

app = FastAPI(lifespan=lifespan)

app.mount("/static", StaticFiles(directory="/static"), name="static")

app.exception_handler(Exception)(handlers.generic_exception_handler)
app.exception_handler(NotAuthenticatedException)(handlers.not_authenticated_handler)

app.add_middleware(context.ContextMiddleware)
app.middleware("http")(mw.timing_middleware)

app.include_router(routes.pages.user.router)
app.include_router(routes.pages.home.router)
app.include_router(routes.pages.auth.router)
app.include_router(routes.pages.files.router)

app.include_router(routes.api.base.router)
app.include_router(routes.api.websockets.router)

app.include_router(routes.resources.router)

app.include_router(routes.htmx.auth_htmx.router)
app.include_router(routes.htmx.search_htmx.router)