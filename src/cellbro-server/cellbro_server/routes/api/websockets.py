import asyncio
from fastapi import WebSocket, APIRouter, WebSocketDisconnect

from ...core.cache import worker_redis
from ...core.templates import templates


router = APIRouter(prefix="/api/ws", tags=["websockets", "api"])

@router.websocket("/output/{task_id}")
async def subscribe_to_worker_messages(websocket: WebSocket, task_id: str):
    await websocket.accept()
    pubsub = worker_redis.client.pubsub()
    output_channel = f"task:{task_id}"
    await pubsub.subscribe(output_channel, "current_task", "step_completed", "task_completed", "task_started", "event_triggered")

    async def listen_redis():
        while True:
            message = await pubsub.get_message(ignore_subscribe_messages=True, timeout=0.1)
            if message and message["type"] == "message":
                if message["channel"] == "current_task":
                    current_task = message["data"]
                    status_output = templates.get_template("components/mini/worker-status.html").render(current_task=current_task)
                    await websocket.send_text(status_output)
                elif message["channel"] == output_channel:
                    output = templates.get_template("components/stdout.html").render(messages=[message["data"]])
                    await websocket.send_text(output)
                # elif message["channel"] == "step_completed":
                #     step_completed = message["data"]
                #     messages = [{"message": f"Step {step_completed} completed!", "category": "success"}]
                #     print(messages, flush=True)
                #     notification = templates.get_template("components/mini/notification.html").render(messages=messages)
                #     await websocket.send_text(notification)
                elif message["channel"] == "task_completed":
                    task_completed = message["data"]
                    messages = [{"message": f"Task {task_completed} completed!", "category": "success"}]
                    notification = templates.get_template("components/mini/notification.html").render(messages=messages)
                    await websocket.send_text(notification)
                elif message["channel"] == "task_started":
                    task_started = message["data"]
                    messages = [{"message": f"Task {task_started} started!", "category": "info"}]
                    notification = templates.get_template("components/mini/notification.html").render(messages=messages)
                    await websocket.send_text(notification) 
                elif message["channel"] == "event_triggered":
                    event_triggered = message["data"]
                    print(f"Event triggered: {event_triggered}", flush=True)
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