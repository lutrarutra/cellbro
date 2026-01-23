import os
from cellbro_db import DBHandler, types

from .. import celery_app, tools


def connect() -> DBHandler:
    db = DBHandler(auto_commit=True)
    db.connect(
        user=os.environ["POSTGRES_USER"],
        password=os.environ["POSTGRES_PASSWORD"],
        host="postgres",
        port=os.environ["POSTGRES_PORT"],
        db=os.environ["POSTGRES_DB"],
    )
    return db


@tools.wrapper.worker_task("io.read_h5ad", complete_steps=types.ChecklistStep.LOAD, trigger_events="dataset-updated")
def read_h5ad(self, file_path: str):
    db = connect()
    from . import io
    with db as session:
        io.read_h5ad(db=session, file_path=file_path)


@tools.wrapper.worker_task("qc", complete_steps=types.ChecklistStep.QC, trigger_events="dataset-updated")
def qc(self, mt_prefix: str, ribo_prefixes: list[str], hb_pattern: str, percent_top: int):
    db = connect()
    import scanpy as sc

    with db as session:
        celery_app.adata.var["mt"] = celery_app.adata.var_names.str.startswith(mt_prefix)
        celery_app.adata.var["ribo"] = celery_app.adata.var_names.str.startswith(tuple(ribo_prefixes))
        celery_app.adata.var["hb"] = celery_app.adata.var_names.str.contains(hb_pattern, regex=True)

        sc.pp.calculate_qc_metrics(
            celery_app.adata,
            qc_vars=["mt", "ribo", "hb"],
            percent_top=[percent_top],
            inplace=True,
            log1p=True
        )

        tools.dataset.reset_dataset(session, celery_app.adata)


