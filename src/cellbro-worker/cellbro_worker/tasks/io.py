from cellbro_db import SyncSession

from .. import celery_app, tools
from ..tools import worker_redis


def read_h5ad(db: SyncSession, file_path: str):
    import anndata as ad
    celery_app.adata = ad.read_h5ad(file_path)
    celery_app.adata.obs_names_make_unique()
    celery_app.adata.var_names_make_unique()

    for key in worker_redis.scan_iter("step:*"):
        worker_redis.delete(key)

    tools.dataset.reset_dataset(db, celery_app.adata)
