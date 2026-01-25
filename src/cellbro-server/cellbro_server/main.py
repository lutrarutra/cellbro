from fastapi import FastAPI
from fastapi.staticfiles import StaticFiles

from .core import middleware as mw, context
from .core.lifespan import lifespan
from .core import exceptions as exc
from . import routes

app = FastAPI(lifespan=lifespan)

app.mount("/static", StaticFiles(directory="/static"), name="static")

app.exception_handler(Exception)(exc.generic_exception_handler)
app.exception_handler(exc.NotAuthenticatedException)(exc.NotAuthenticatedException.handler)
app.exception_handler(exc.DatasetNotLoadedException)(exc.DatasetNotLoadedException.handler)
app.exception_handler(exc.QCNotCompletedException)(exc.QCNotCompletedException.handler)


app.add_middleware(context.ContextMiddleware)
app.middleware("http")(mw.timing_middleware)

app.include_router(routes.pages.user.router)
app.include_router(routes.pages.dashboard.router)
app.include_router(routes.pages.auth.router)
app.include_router(routes.pages.files.router)
app.include_router(routes.pages.qc.router)

app.include_router(routes.api.base.router)
app.include_router(routes.api.websockets.router)

app.include_router(routes.resources.router)

app.include_router(routes.htmx.auth_htmx.router)
app.include_router(routes.htmx.search_htmx.router)
app.include_router(routes.htmx.data_htmx.router)
app.include_router(routes.htmx.steps_htmx.router)
app.include_router(routes.htmx.plots_htmx.router)