from cellbro_db import types, SyncSession

from taskiq import TaskiqDepends

from .. import tools, db_session
from ..tools.wrappers import worker_task
from ..broker import task_broker

@task_broker.task
@worker_task(complete_steps=types.ChecklistStep.QC)
async def run(mt_prefix: str, ribo_prefixes: list[str], hb_pattern: str, percent_top: int, session: SyncSession = TaskiqDepends(db_session)):
    import scanpy as sc
    import anndata as ad

    adata = ad.read_zarr("/app/data/dataset.zarr")
    adata.var["mt"] = adata.var_names.str.startswith(mt_prefix)
    adata.var["ribo"] = adata.var_names.str.startswith(tuple(ribo_prefixes))
    adata.var["hb"] = adata.var_names.str.contains(hb_pattern, regex=True)

    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt", "ribo", "hb"],
        percent_top=[percent_top],
        inplace=True,
        log1p=True
    )

    import zarr
    store = zarr.open("/app/data/dataset.zarr", mode="r+")
    ad.io.write_elem(store, "obs", adata.obs)
    ad.io.write_elem(store, "var", adata.var)

    tools.dataset.reset_dataset(session, adata)

@task_broker.task
@worker_task()
async def plot_total_counts_histogram(plot_id: str):
    from bokeh.plotting import figure
    from bokeh.themes import built_in_themes
    from bokeh.document import Document
    from bokeh.core import types
    from bokeh.embed import json_item
    import numpy as np
    import anndata as ad

    adata = ad.read_zarr("/app/data/dataset.zarr")
    counts = adata.obs["total_counts"].var
    hist, edges = np.histogram(counts, bins=50)  # type: ignore

    p = figure(
        x_axis_label="Total Counts", 
        y_axis_label="Number of Cells", 
        sizing_mode="stretch_both", 
        output_backend="webgl"
    )

    p.quad(
        top=hist,
        bottom=0,
        left=edges[:-1],
        right=edges[1:],
        fill_color="navy",
        line_color="white", # Typo fixed
        alpha=0.5,
    )

    doc = Document()
    doc.add_root(p)
    doc.theme = built_in_themes['dark_minimal']

    return json_item(p, types.ID(plot_id))