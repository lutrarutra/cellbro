from .FlashCache import FlashCache
from .SessionCache import SessionCache
from .RedisManager import RedisManager
from .config import settings

session_cache = SessionCache(expiration=settings.SESSION_EXPIRE_SECONDS)
flash_cache = FlashCache()
worker_output_cache = RedisManager()