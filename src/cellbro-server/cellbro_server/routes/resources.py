from pathlib import Path

from fastapi import APIRouter
from fastapi.responses import FileResponse

from ..core import exceptions as exc

router = APIRouter(tags=["resources"])


# @router.get("/static/{resource_path:path}")
# async def get_static_resource(resource_path: Path):
#     file_path = Path("/static/") / resource_path
#     if not file_path.exists() or not file_path.is_file():
#         raise exc.ResourceNotFound(file_path.as_posix())

#     return FileResponse(file_path)