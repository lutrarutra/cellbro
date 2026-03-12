import sqlalchemy as sa
import pandas as pd
import anndata as ad

from cellbro_db import SyncSession, models, types, queries

def clean_db(db: SyncSession):
    db.execute(sa.delete(models.Observation))
    db.execute(sa.delete(models.Feature))
    db.execute(sa.delete(models.Layer))
    db.execute(sa.delete(models.Variable))


def reset_dataset(db: SyncSession, adata: ad.AnnData):
    clean_db(db)
    
    objects_to_add = []
    for idx in adata.obs.index:
        objects_to_add.append(queries.observation.create(name=str(idx)))

    for idx in adata.var.index:
        feature = queries.feature.create(identifier=str(idx), name=str(idx))
        objects_to_add.append(feature)

    objects_to_add.extend([
        queries.layer.create(name=layer_name, dtype=str(adata.layers[layer_name].dtype))
        for layer_name in adata.layers.keys()
    ])
    if adata.X is not None:
        objects_to_add.append(queries.layer.create(name="X", dtype=str(adata.X.dtype)))

    for var_name in adata.var.columns:
        dtype = adata.var[var_name].dtype

        if pd.api.types.is_integer_dtype(dtype):
            var_type = types.VariableType.INTEGER
        elif pd.api.types.is_float_dtype(dtype):
            var_type = types.VariableType.FLOAT
        elif pd.api.types.is_bool_dtype(dtype):
            var_type = types.VariableType.BOOLEAN
        elif isinstance(dtype, pd.CategoricalDtype):
            var_type = types.VariableType.CATEGORICAL
        elif pd.api.types.is_string_dtype(dtype):
            var_type = types.VariableType.STRING
        else:
            var_type = types.VariableType.STRING

        objects_to_add.append(queries.variable.create(
            name=var_name,
            layer=types.AnnDataLayerType.VAR,
            type=var_type,
        ))

    for var_name in adata.obs.columns:
        dtype = adata.obs[var_name].dtype

        if pd.api.types.is_integer_dtype(dtype):
            var_type = types.VariableType.INTEGER
        elif pd.api.types.is_float_dtype(dtype):
            var_type = types.VariableType.FLOAT
        elif pd.api.types.is_bool_dtype(dtype):
            var_type = types.VariableType.BOOLEAN
        elif isinstance(dtype, pd.CategoricalDtype):
            var_type = types.VariableType.CATEGORICAL
        elif pd.api.types.is_string_dtype(dtype):
            var_type = types.VariableType.STRING
        else:
            var_type = types.VariableType.STRING

        objects_to_add.append(queries.variable.create(
            name=var_name,
            layer=types.AnnDataLayerType.OBS,
            type=var_type,
        ))

    db.add_all(objects_to_add)
    db.commit()