from .FlashCache import FlashCache
from .SessionCache import SessionCache
from .config import settings

session_cache = SessionCache(expiration=settings.SESSION_EXPIRE_SECONDS)
flash_cache = FlashCache()