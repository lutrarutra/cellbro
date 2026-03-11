from redis.asyncio import Redis

class SessionCache:
    expiration: int | None = None
    client: Redis = None # type: ignore

    def __init__(self, expiration: int | None = None):
        super().__init__()
        self.expiration = expiration

    def connect(self, host: str, port: int, db: int, decode_responses: bool = True):
        print(f"Connecting to Redis at {host}:{port}, db={db}")
        self.client = Redis(host=host, port=port, db=db, decode_responses=decode_responses)

    async def close(self):
        if self.client:
            await self.client.close()
            self.client = None  # type: ignore

    def _redis_key(self, sid: str) -> str:
        return sid

    async def set(self, sid: str, used_id: int) -> None:
        await self.client.set(self._redis_key(sid), str(used_id), ex=self.expiration)
        
    async def get(self, sid: str) -> int | None:
        if (res := await self.client.get(self._redis_key(sid))) is not None:
            return int(res)
        return None

    async def delete(self, sid: str) -> None:
        await self.client.delete(self._redis_key(sid))

    async def delete_user_sessions(self, user_id: int) -> int:
        lua_script = """
        local keys_to_delete = {}
        -- Warning: KEYS * is slow on large databases
        local session_keys = redis.call('KEYS', '*') 
        local target_id = tostring(ARGV[1])
        
        for i, key in ipairs(session_keys) do
            local session_val = redis.call('GET', key)
            if session_val == target_id then
                table.insert(keys_to_delete, key)
            end
        end
        
        if #keys_to_delete > 0 then
            return redis.call('DEL', unpack(keys_to_delete))
        end
        return 0
""" 
        return int(await self.client.eval(lua_script, 0, user_id))  # type: ignore