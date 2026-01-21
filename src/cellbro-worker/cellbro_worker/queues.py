from celery.result import AsyncResult

def read_h5ad(file_path: str) -> AsyncResult:
    from .tasks import read_h5ad
    return read_h5ad.apply_async(args=[file_path])  # type: ignore