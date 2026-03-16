from cellbro_db import types, SyncSession

from taskiq import TaskiqDepends

from .. import tools, db_session
from ..broker import task_broker

@task_broker.task
@tools.wrappers.worker_task(complete_steps=types.ChecklistStep.LOAD, trigger_events="dataset_loaded")
async def read_h5ad(file_path: str, session: SyncSession = TaskiqDepends(db_session)):
    import anndata as ad
    ad.settings.allow_write_nullable_strings = True

    adata = ad.read_h5ad(file_path, backed="r")
    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    adata.write_zarr("/app/data/dataset.zarr")
    
    tools.dataset.reset_dataset(session, adata)