import os
import json

from cellbro_db import DBHandler, types

from .. import CellBroWorker
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


@tools.wrapper.worker_task("io.read_h5ad", complete_steps=types.ChecklistStep.LOAD, trigger_events="dataset-updated", notify=True, read_resources=[], write_resources=["X", "var", "obs", "uns", "layers", "obsm", "varm"])
def read_h5ad(self, file_path: str):
    print("##$$$$")
    db = connect()
    from . import io
    with db as session:
        io.read_h5ad(db=session, file_path=file_path)
    print(celery_app.adata)
    print(self.app.adata, flush=True)



@tools.wrapper.worker_task("qc", complete_steps=types.ChecklistStep.QC, notify=True, trigger_events="dataset-updated", read_resources=["X"], write_resources=["var", "obs"])
def qc(self, mt_prefix: str, ribo_prefixes: list[str], hb_pattern: str, percent_top: int):
    db = connect()
    import scanpy as sc

    self: CellBroWorker = self.app

    with db as session:
        self.adata.var["mt"] = self.adata.var_names.str.startswith(mt_prefix)
        self.adata.var["ribo"] = self.adata.var_names.str.startswith(tuple(ribo_prefixes))
        self.adata.var["hb"] = self.adata.var_names.str.contains(hb_pattern, regex=True)

        sc.pp.calculate_qc_metrics(
            self.adata,
            qc_vars=["mt", "ribo", "hb"],
            percent_top=[percent_top],
            inplace=True,
            log1p=True
        )

        tools.dataset.reset_dataset(session, self.adata)


@tools.wrapper.worker_task("plot.test_plot", notify=False, read_resources=["X"])
def test_plot(self, plot_id: str):
    from bokeh.plotting import figure
    from bokeh.themes import built_in_themes
    from bokeh.document import Document
    from bokeh.embed import json_item

    import time
    time.sleep(1)

    p = figure(x_axis_label="X", y_axis_label="Y", sizing_mode="stretch_both", match_aspect=True, output_backend="webgl")
    p.line([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], line_width=2, color="navy")

    doc = Document()
    doc.add_root(p)
    doc.theme = built_in_themes['dark_minimal']

    tools.worker_redis.publish(f'plot:{plot_id}', json.dumps(json_item(p, plot_id)))  # type: ignore


@tools.wrapper.worker_task("plot.total_counts_histogram", read_resources=["X", "obs", "var"])
def plot_total_counts_histogram(self, plot_id: str):
    from bokeh.plotting import figure
    from bokeh.themes import built_in_themes
    from bokeh.document import Document
    from bokeh.embed import json_item
    import anndata as ad
    import numpy as np

    adata: ad.AnnData = self.app.adata

    counts = adata.obs["total_counts"].values
    hist, edges = np.histogram(counts, bins=50)

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

    tools.worker_redis.publish(f'plot:{plot_id}', json.dumps(json_item(p, plot_id)))  # type: ignore