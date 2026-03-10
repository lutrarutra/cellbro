from cellbro_db import types, SyncSession

from taskiq import TaskiqDepends
from .. import tools, broker, db_session


@tools.wrapper.worker_task(
    task_name="io:read_h5ad",
    read_resources=["dataset"],
    write_resources=["dataset"],
    complete_steps=types.ChecklistStep.LOAD,
    notify=True
)
def read_h5ad(file_path: str, session: SyncSession = TaskiqDepends(db_session)):
    import anndata as ad

    adata = ad.read_h5ad(file_path)
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    
    tools.dataset.reset_dataset(session, adata)