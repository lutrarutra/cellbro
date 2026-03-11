import json
from typing import Literal
from cellbro_db import types
from redis.asyncio import Redis as AsyncRedis
from redis import Redis as SyncRedis

class AsyncMessageRedis:
    client: AsyncRedis = None # type: ignore
    
    def connect(self, host: str, port: int, db: int, decode_responses: bool = True):
        print(f"[Async] Connecting to Redis at {host}:{port}, db={db}")
        self.client = AsyncRedis(host=host, port=port, db=db, decode_responses=decode_responses)
        
    async def log(self, task_id: str, text: str, category: str = "log"):
        await self.client.publish("stdout", json.dumps({"task_id": task_id, "text": text, "category": category}))
        
    async def error(self, task_id: str, text: str):
        await self.log(task_id, text, category="error")
        
    async def warning(self, task_id: str, text: str):
        await self.log(task_id, text, category="warning")

    async def notify(self, text: str, category: Literal["log", "error", "warning"] = "log"):
        await self.client.publish("notify", json.dumps({"text": text, "category": category}))

    async def complete_task(self, task_id: str):
        await self.client.publish(f"task:{task_id}", "completed")
        
    async def step(self, step_name: str | int, status: str):
        await self.client.set(f"step:{step_name}", status)
        await self.client.publish("step_completed", step_name)
        
    async def event(self, event_name: str):
        await self.client.publish("event_triggered", event_name)
    
    async def close(self):
        if self.client:
            await self.client.aclose()
            self.client = None  # type: ignore
            
    async def get_completed_steps(self) -> list[types.ChecklistStep]:
        completed_steps = []
        for step in types.ChecklistStep:
            if (await self.client.get(f"step:{step}") == "completed"):
                completed_steps.append(step)
        return completed_steps
    
    async def is_step_completed(self, step: types.ChecklistStep) -> bool:
        return await self.client.get(f"step:{step}") == "completed"
    
    async def get_current_task(self) -> str | None:
        return await self.client.get("current_task")
    
    async def get_number_of_running_tasks(self) -> int:
        count = await self.client.get("running_tasks")
        return int(count) if count is not None else 0


class SyncMessageRedis:
    client: SyncRedis = None # type: ignore
    
    def connect(self, host: str, port: int, db: int, decode_responses: bool = True):
        print(f"[Sync] Connecting to Redis at {host}:{port}, db={db}")
        self.client = SyncRedis(host=host, port=port, db=db, decode_responses=decode_responses)
        
    def log(self, task_id: str, text: str, category: str = "log"):
        self.client.publish("stdout", json.dumps({"task_id": task_id, "text": text, "category": category}))
        
    def error(self, task_id: str, text: str):
        self.log(task_id, text, category="error")
        
    def warning(self, task_id: str, text: str):
        self.log(task_id, text, category="warning")

    def notify(self, text: str, category: Literal["log", "error", "warning"] = "log"):
        self.client.publish("notify", json.dumps({"text": text, "category": category}))
        
    def complete_task(self, task_id: str):
        self.client.publish(f"task:{task_id}", "completed")
        
    def step(self, step_name: str | int, status: str):
        self.client.set(f"step:{step_name}", status)
        self.client.publish("step_completed", step_name)
        
    def event(self, event_name: str):
        self.client.publish("event_triggered", event_name)
    
    def close(self):
        if self.client:
            self.client.close()
            self.client = None  # type: ignore
            
    def get_completed_steps(self) -> list[types.ChecklistStep]:
        completed_steps = []
        for step in types.ChecklistStep:
            if (self.client.get(f"step:{step}") == "completed"):
                completed_steps.append(step)
        return completed_steps
    
    def is_step_completed(self, step: types.ChecklistStep) -> bool:
        return self.client.get(f"step:{step}") == "completed"
    
    def get_current_task(self) -> str | None:
        return self.client.get("current_task")  # type: ignore
    
    def get_number_of_running_tasks(self) -> int:
        count = self.client.get("running_tasks")
        return int(count) if count is not None else 0  # type: ignore