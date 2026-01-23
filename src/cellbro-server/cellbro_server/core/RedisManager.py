from cellbro_db import models, types

from redis.asyncio import Redis

class RedisManager:
    client: Redis = None # type: ignore
    
    async def connect(self, host: str, port: int, db: int, decode_responses: bool = True):
        print(f"Connecting to Redis at {host}:{port}, db={db}")
        self.client = Redis(host=host, port=port, db=db, decode_responses=decode_responses)
    
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
        status = await self.client.get(f"step:{step}")
        return status == "completed"
    
    async def get_current_task(self) -> str | None:
        return await self.client.get("current_task")