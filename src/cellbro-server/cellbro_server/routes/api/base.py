from fastapi import APIRouter, Depends
from cellbro_worker import tasks

router = APIRouter(tags=["base", "api"])

from ...core.cache import flash_cache, worker_redis
from ...core.context import ctx
from ...core.responses import html_response
from ...core.dependencies import get_sid


@router.get("/status")
async def status():
    return {"status": "ok"}


@router.get("/worker_status")
async def worker_status():
    running_task_count = await worker_redis.get_number_of_running_tasks()
    return await html_response("components/mini/worker-status.html", running_task_count=running_task_count)

@router.get("/flash-messages")
async def retrieve_flash_messages(sid: str = Depends(get_sid)):
    if (messages := await flash_cache.consume(sid)):
        return {"messages": messages}
    return {"messages": []}


@router.post("/read_data")
async def read_data():
    task_id = await tasks.io.read_h5ad.kiq("/app/data/pbmc3k.h5ad")
    return await html_response("components/mini/loading.html", task_id=task_id, task_name="Read Data")

@router.get("/timeline")
async def timeline():
    completed_steps = await worker_redis.get_completed_steps()
    return await html_response("timeline.html", completed_steps=completed_steps)