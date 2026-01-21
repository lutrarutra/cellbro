from cellbro_db import SyncSession, models

from .. import celery_app


def read_h5ad(db: SyncSession, file_path: str):
    import anndata as ad
    celery_app.adata = ad.read_h5ad(file_path)
    celery_app.adata.obs_names_make_unique()
    celery_app.adata.var_names_make_unique()
    observations = []

    for idx in celery_app.adata.obs.index:
        observations.append(models.Observation.Create(name=str(idx)))


    db.add_all(observations)
    print(f"Added {len(observations)} observations to the database.")
    del observations

    features = []
    for idx in celery_app.adata.var.index:
        feature = models.Feature.Create(identifier=str(idx), name=str(idx))
        features.append(feature)

    db.add_all(features)
    print(f"Added {len(features)} features to the database.")
    del features
