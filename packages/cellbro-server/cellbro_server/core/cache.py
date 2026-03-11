from .FlashCache import FlashCache
from .SessionCache import SessionCache
from .config import settings
from cellbro_worker.tools.MessageRedis import AsyncMessageRedis

session_cache = SessionCache(expiration=settings.SESSION_EXPIRE_SECONDS)
flash_cache = FlashCache()
worker_redis = AsyncMessageRedis()