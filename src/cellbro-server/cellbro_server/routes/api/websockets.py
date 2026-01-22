import asyncio
from fastapi import WebSocket, APIRouter, WebSocketDisconnect

from ...core.cache import worker_output_cache
from ...core.templates import templates


router = APIRouter(prefix="/api/ws", tags=["websockets", "api"])

@router.websocket("/output/{task_id}")
async def stream_task_stdout(websocket: WebSocket, task_id: str):
    await websocket.accept()
    pubsub = worker_output_cache.client.pubsub()
    await pubsub.subscribe(f"task:{task_id}")

    async def listen_redis():
        while True:
            message = await pubsub.get_message(ignore_subscribe_messages=True, timeout=1.0)
            if message and message["type"] == "message":
                output = templates.get_template("components/stdout.html").render(messages=[message["data"]])
                print(output, flush=True)
                await websocket.send_text(output)

    async def receive_from_client():
        while True:
            await websocket.receive_text()

    listener_task = asyncio.create_task(listen_redis())
    receiver_task = asyncio.create_task(receive_from_client())

    try:
        done, pending = await asyncio.wait(
            {listener_task, receiver_task},
            return_when=asyncio.FIRST_COMPLETED
        )
        
        if listener_task in done:
            listener_task.result()

    except WebSocketDisconnect:
        print("Client disconnected normally")
    except Exception as e:
        print(f"WebSocket error: {e}")
    finally:
        print("Cleaning up...")
        listener_task.cancel()
        receiver_task.cancel()
        await pubsub.close()