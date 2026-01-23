from celery.result import AsyncResult

from cellbro_db import types

def read_h5ad(file_path: str) -> AsyncResult:
    from .tasks import read_h5ad
    return read_h5ad.apply_async(args=[file_path])  # type: ignore

def qc(mt_prefix: str, ribo_prefixes: list[str], hb_pattern: str, percent_top: int) -> AsyncResult:
    from .tasks import qc
    return qc.apply_async(args=[mt_prefix, ribo_prefixes, hb_pattern, percent_top])  # type: ignore