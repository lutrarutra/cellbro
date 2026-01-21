from fastapi import APIRouter, Depends

router = APIRouter(tags=["base", "api"])

from ...core.cache import flash_cache
from ...core.dependencies import get_sid


@router.get("/status")
async def status():
    return {"status": "ok"}

@router.get("/flash-messages")
async def retrieve_flash_messages(sid: str = Depends(get_sid)):
    if (messages := await flash_cache.consume(sid)):
        return {"messages": messages}
    return {"messages": []}