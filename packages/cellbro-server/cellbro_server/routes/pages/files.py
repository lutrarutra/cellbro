import os
from pathlib import Path
from fastapi import Depends, APIRouter
from fastapi.concurrency import run_in_threadpool

from cellbro_db.core.session import AsyncSession
from cellbro_db import models

from ...core import exceptions as exc
from ...core.dependencies import db_session, get_user
from ...core.responses import html_response

router = APIRouter(prefix="/files", tags=["files", "view"])


class PathEntry:
    def __init__(self, name: str, path: str, is_dir: bool, size: int | None = None, modified: float | None = None):
        self.name = name
        self.path = path
        self.is_dir = is_dir
        self.size = size
        self.modified = modified

def get_directory_contents_sync(target_path: str):
    entries: list[PathEntry] = []
    
    try:
        with os.scandir(target_path) as it:
            for entry in it:
                stat = entry.stat()
                entries.append(PathEntry(
                    name=entry.name,
                    path=entry.path,
                    is_dir=entry.is_dir(),
                    size=stat.st_size,
                    modified=stat.st_mtime
                ))
        return sorted(entries, key=lambda e: e.name)
    except PermissionError:
        return {"error": "Permission Denied"}


@router.get("/")
async def browse(path: str | None = None, user: models.User = Depends(get_user), db: AsyncSession = Depends(db_session)):
    if not path:
        path = os.path.abspath(os.sep)
    paths = await run_in_threadpool(get_directory_contents_sync, path)
    parent_dir = str(Path(path).parent) if Path(path).parent != Path(path) else None
    return await html_response("views/files.html", paths=paths, current_path=path, parent_dir=parent_dir)