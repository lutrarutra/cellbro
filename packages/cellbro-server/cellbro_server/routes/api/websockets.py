import asyncio
import json
from fastapi import WebSocket, APIRouter, WebSocketDisconnect

from ...core.cache import worker_redis
from ...core.templates import templates


router = APIRouter(prefix="/api/ws", tags=["websockets", "api"])

@router.websocket("/worker-messages")
async def subscribe_to_worker_messages(websocket: WebSocket):
    await websocket.accept()
    pubsub = worker_redis.client.pubsub()
    await pubsub.psubscribe("task_status:*")
    await pubsub.subscribe("stdout", "step_completed", "task_completed", "task_started", "event_triggered", "notify")

    async def listen_redis():
        while True:
            message = await pubsub.get_message(ignore_subscribe_messages=True, timeout=0.1)  # type: ignore
            if message and message["type"] in ["message", "pmessage"]:
                if message["channel"].startswith("task_status:"):
                    running_task_count = await worker_redis.get_number_of_running_tasks()
                    status_output = templates.get_template("components/mini/worker-status.html").render(running_task_count=running_task_count)
                    await websocket.send_text(status_output)
                elif message["channel"] == "stdout":
                    output = templates.get_template("components/mini/stdout.html").render(message=json.loads(message["data"]))
                    await websocket.send_text(output)
                elif message["channel"] == "notify":
                    notification = templates.get_template("components/mini/notification.html").render(notification_data = message["data"])
                    await websocket.send_text(notification)
                elif message["channel"] == "event_triggered":
                    event_triggered = message["data"]
                    notification = templates.get_template("components/mini/htmx-trigger.html").render(event=event_triggered)
                    await websocket.send_text(notification)

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


@router.websocket("/loading/{task_id}")
async def loading_listener(websocket: WebSocket, task_id: str):
    await websocket.accept()
    pubsub = worker_redis.client.pubsub()
    await pubsub.subscribe(f"task:{task_id}")

    async def listen():
        while True:
            message = await pubsub.get_message(ignore_subscribe_messages=True, timeout=0.1)  # type: ignore
            if message and message["type"] == "message":
                print(f"Received loading message on channel {message['channel']}: {message['data']}", flush=True)
                if message["channel"] == f"task:{task_id}":
                    redirect_to = await worker_redis.client.get(f"task_redirect:{task_id}")
                    if (status := message["data"]) == "completed":
                        if redirect_to:
                            response = f'<div id="global-loader" hx-swap-oob="true" hx-get="{redirect_to}" hx-trigger="load" hx-target="#content-container"></div>'
                        else:
                            response = templates.get_template("components/mini/loading.html").render(done=True)
                    elif status == "failed":
                        if redirect_to:
                            response = f'<div id="global-loader" hx-swap-oob="true" hx-get="{redirect_to}" hx-trigger="load" hx-target="#content-container"></div>'
                        else:
                            response = templates.get_template("components/mini/loading.html").render(done=True)
                    else:
                        continue
                    await websocket.send_text(response)
                    break

            await asyncio.sleep(0.1)

    async def receive_from_client():
        while True:
            await websocket.receive_text()

    listener_task = asyncio.create_task(listen())
    receiver_task = asyncio.create_task(receive_from_client())

    try:
        done, pending = await asyncio.wait({listener_task, receiver_task}, return_when=asyncio.FIRST_COMPLETED)
        
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