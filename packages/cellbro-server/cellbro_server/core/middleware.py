import time

from typing import Callable, Awaitable
from fastapi import Request, Response

async def timing_middleware(request: Request, call_next: Callable[[Request], Awaitable[Response]]):
    start_time = time.time()
    response = await call_next(request)

    process_time = time.time() - start_time
    response.headers["X-Process-Time"] = str(process_time)
    return response