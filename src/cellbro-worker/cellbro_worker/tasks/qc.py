# from cellbro_db import types

# from .. import tools, CellBroTask


# @tools.wrapper.worker_task("qc.run", complete_steps=types.ChecklistStep.QC, notify=True, trigger_events="dataset-updated", read_resources=["X"], write_resources=["var", "obs"])
# def run(self: CellBroTask, mt_prefix: str, ribo_prefixes: list[str], hb_pattern: str, percent_top: int):
#     import scanpy as sc

#     import time
#     time.sleep(5)

#     with self.db as session:
#         self.adata.var["mt"] = self.adata.var_names.str.startswith(mt_prefix)
#         self.adata.var["ribo"] = self.adata.var_names.str.startswith(tuple(ribo_prefixes))
#         self.adata.var["hb"] = self.adata.var_names.str.contains(hb_pattern, regex=True)

#         sc.pp.calculate_qc_metrics(
#             self.adata,
#             qc_vars=["mt", "ribo", "hb"],
#             percent_top=[percent_top],
#             inplace=True,
#             log1p=True
#         )

#         tools.dataset.reset_dataset(session, self.adata)

# @tools.wrapper.worker_task("qc.plot.total_counts_histogram", read_resources=["X", "obs", "var"])
# def plot_total_counts_histogram(self: CellBroTask, plot_id: str):
#     from bokeh.plotting import figure
#     from bokeh.themes import built_in_themes
#     from bokeh.document import Document
#     from bokeh.core import types
#     from bokeh.embed import json_item
#     import numpy as np
#     import json

#     counts = self.adata.obs["total_counts"].values
#     hist, edges = np.histogram(counts, bins=50)  # type: ignore

#     p = figure(
#         x_axis_label="Total Counts", 
#         y_axis_label="Number of Cells", 
#         sizing_mode="stretch_both", 
#         output_backend="webgl"
#     )

#     p.quad(
#         top=hist,
#         bottom=0,
#         left=edges[:-1],
#         right=edges[1:],
#         fill_color="navy",
#         line_color="white", # Typo fixed
#         alpha=0.5,
#     )

#     doc = Document()
#     doc.add_root(p)
#     doc.theme = built_in_themes['dark_minimal']

#     self.r.publish(f'plot:{plot_id}', json.dumps(json_item(p, types.ID(plot_id))))