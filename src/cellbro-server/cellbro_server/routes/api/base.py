from multiprocessing import queues
from fastapi import APIRouter, Depends
from cellbro_worker import queues

router = APIRouter(tags=["base", "api"])

from ...core.cache import flash_cache, worker_output_cache
from ...core.responses import html_response
from ...core.dependencies import get_sid


@router.get("/status")
async def status():
    return {"status": "ok"}


@router.get("/worker_status")
async def worker_status():
    status = await worker_output_cache.client.get("status")
    print(f"Worker status: {status}", flush=True)
    busy = status == "busy"
    return await html_response("components/worker-status.html", busy=busy)

@router.get("/flash-messages")
async def retrieve_flash_messages(sid: str = Depends(get_sid)):
    if (messages := await flash_cache.consume(sid)):
        return {"messages": messages}
    return {"messages": []}


@router.post("/read_data")
async def read_data():
    task_id = queues.read_h5ad("/app/data/pbmc3k_raw.h5ad")
    return await html_response(status=204)