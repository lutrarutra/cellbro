from redis.asyncio import Redis

class FlashCache:
    PRIORITY = ["error", "warning", "info", "success"]
    client: Redis
    
    def __init__(self, expiration: int = 3600):
        self.expiration = expiration

    async def connect(self, host: str, port: int, db: int):
        self.client = Redis(host=host, port=port, db=db, decode_responses=True)

    def _redis_key(self, sid: str, category: str) -> str:
        return f"flash:{sid}:{category}"
    
    async def add(self, sid: str, category: str, message: str) -> None:
        key = self._redis_key(sid, category)
        async with self.client.pipeline(transaction=True) as pipe:
            await pipe.rpush(key, message)
            await pipe.expire(key, self.expiration)
            await pipe.execute()

    async def consume(self, sid: str) -> list[tuple[str, str]]:
        res = []
        for category in self.PRIORITY:
            key = self._redis_key(sid, category)
            messages: list[str] = await self.client.lrange(key, 0, -1)
            if messages:
                await self.client.delete(key)
                res.extend((category, msg) for msg in messages)
        return res
            
    async def close(self) -> None:
        await self.client.close()