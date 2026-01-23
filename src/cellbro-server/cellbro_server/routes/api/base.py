from multiprocessing import queues
from fastapi import APIRouter, Depends
from cellbro_worker import queues

router = APIRouter(tags=["base", "api"])

from ...core.cache import flash_cache, worker_redis
from ...core.responses import html_response
from ...core.dependencies import get_sid


@router.get("/status")
async def status():
    return {"status": "ok"}


@router.get("/worker_status")
async def worker_status():
    current_task = await worker_redis.client.get("current_task")
    return await html_response("components/mini/worker-status.html", current_task=current_task)

@router.get("/flash-messages")
async def retrieve_flash_messages(sid: str = Depends(get_sid)):
    if (messages := await flash_cache.consume(sid)):
        return {"messages": messages}
    return {"messages": []}


@router.post("/read_data")
async def read_data():
    task_id = queues.read_h5ad("/app/data/pbmc3k.h5ad")
    return await html_response(status=204)

@router.get("/timeline")
async def timeline():
    completed_steps = await worker_redis.get_completed_steps()
    return await html_response("timeline.html", completed_steps=completed_steps)