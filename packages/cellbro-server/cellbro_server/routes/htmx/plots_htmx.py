from fastapi import Depends, APIRouter, HTTPException, Response
import sqlalchemy as sa
import asyncio
import json

router = APIRouter(prefix="/htmx/plots", tags=["plots", "htmx"])

from cellbro_db import models
from cellbro_db.core.session import AsyncSession
from cellbro_worker import queues

from ...core.context import ctx
from ...core import responses
from ... import forms, logic
from ...core.cache import worker_redis


@router.get("/plot")
async def get_plot(plot_id: str) -> Response:
    pubsub = worker_redis.client.pubsub()
    await pubsub.subscribe(f"plot:{plot_id}")

    queues.test_plot(plot_id)

    async def wait_for_response():
        async for message in pubsub.listen():
            if message['type'] == 'message' and message['channel'] == f'plot:{plot_id}':
                return message['data']
    
    try:
        return Response(await asyncio.wait_for(wait_for_response(), timeout=10.0), media_type="application/json")
    except asyncio.TimeoutError:
        raise HTTPException(status_code=504, detail="Worker timed out")
    finally:
        await pubsub.unsubscribe(f"plot:{plot_id}")
        await pubsub.close()


@router.get("/total_counts_histogram")
async def get_total_counts_histogram(plot_id: str) -> Response:
    pubsub = worker_redis.client.pubsub()
    await pubsub.subscribe(f"plot:{plot_id}")

    queues.plot_total_counts_histogram(plot_id)

    async def wait_for_response():
        async for message in pubsub.listen():
            if message['type'] == 'message' and message['channel'] == f'plot:{plot_id}':
                return message['data']
    
    try:
        return Response(await asyncio.wait_for(wait_for_response(), timeout=10.0), media_type="application/json")
    except asyncio.TimeoutError:
        raise HTTPException(status_code=504, detail="Worker timed out")
    finally:
        await pubsub.unsubscribe(f"plot:{plot_id}")
        await pubsub.close()