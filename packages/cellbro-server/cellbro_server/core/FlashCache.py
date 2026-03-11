from redis.asyncio import Redis

class FlashCache():
    client: Redis = None # type: ignore
    PRIORITY = ["error", "warning", "info", "success"]
    
    def __init__(self, expiration: int = 3600):
        super().__init__()
        self.expiration = expiration

    def connect(self, host: str, port: int, db: int, decode_responses: bool = True):
        print(f"Connecting to Redis at {host}:{port}, db={db}")
        self.client = Redis(host=host, port=port, db=db, decode_responses=decode_responses)

    async def close(self):
        if self.client:
            await self.client.aclose()
            self.client = None  # type: ignore

    def _redis_key(self, sid: str, category: str) -> str:
        return f"flash:{sid}:{category}"
    
    async def add(self, sid: str, category: str, message: str) -> None:
        key = self._redis_key(sid, category)
        async with self.client.pipeline(transaction=True) as pipe:
            await pipe.rpush(key, message)  # type: ignore
            await pipe.expire(key, self.expiration)
            await pipe.execute()

    async def consume(self, sid: str) -> list[tuple[str, str]]:
        res = []
        for category in self.PRIORITY:
            key = self._redis_key(sid, category)
            messages: list[str] = await self.client.lrange(key, 0, -1)  # type: ignore
            if messages:
                await self.client.delete(key)
                res.extend((category, msg) for msg in messages)
        return res
