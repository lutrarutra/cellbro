from fastapi import Depends, APIRouter, HTTPException, Response
import sqlalchemy as sa
import asyncio

router = APIRouter(prefix="/htmx/plots", tags=["plots", "htmx"])

from cellbro_worker import tasks

from ...core.context import ctx
from ...core.cache import worker_redis


@router.get("/plot")
async def get_plot(plot_id: str) -> Response:
    task = await tasks.qc.plot_total_counts_histogram.kiq(plot_id)
    pubsub = worker_redis.client.pubsub()
    await pubsub.subscribe(f"result:{task.task_id}")

    async def wait_for_response():
        async for message in pubsub.listen():
            if message['type'] == 'message' and message['channel'] == f'result:{task.task_id}':
                return message['data']
    
    try:
        return Response(await asyncio.wait_for(wait_for_response(), timeout=10.0), media_type="application/json")
    except asyncio.TimeoutError:
        raise HTTPException(status_code=504, detail="Worker timed out")
    finally:
        await pubsub.unsubscribe(f"result:{task.task_id}")
        await pubsub.close()


@router.get("/total_counts_histogram")
async def get_total_counts_histogram(plot_id: str) -> Response:
    task = await tasks.qc.plot_total_counts_histogram.kiq(plot_id)

    pubsub = worker_redis.client.pubsub()
    await pubsub.subscribe(f"result:{task.task_id}")
    async def wait_for_response():
        async for message in pubsub.listen():
            if message['type'] == 'message' and message['channel'] == f'result:{task.task_id}':
                return message['data']
    
    try:
        return Response(await asyncio.wait_for(wait_for_response(), timeout=10.0), media_type="application/json")
    except asyncio.TimeoutError:
        raise HTTPException(status_code=504, detail="Worker timed out")
    finally:
        await pubsub.unsubscribe(f"result:{task.task_id}")
        await pubsub.close()