from cellbro_db import types

from .. import tools, CellBroTask

@tools.wrapper.worker_task("io.read_h5ad", complete_steps=types.ChecklistStep.LOAD, trigger_events="dataset-updated", notify=True, read_resources=[], write_resources=["X", "var", "obs", "uns", "layers", "obsm", "varm"])
def read_h5ad(self: CellBroTask, file_path: str):
    import anndata as ad

    with self.db as session:
        self.adata = ad.read_h5ad(file_path)
        self.adata.obs_names_make_unique()
        self.adata.var_names_make_unique()

        import time
        time.sleep(5)

        for key in self.r.scan_iter("step:*"):
            self.r.delete(key)

        tools.dataset.reset_dataset(session, self.adata)
