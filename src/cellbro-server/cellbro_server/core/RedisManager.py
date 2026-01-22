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