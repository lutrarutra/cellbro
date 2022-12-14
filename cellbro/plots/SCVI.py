import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html
from plotly.subplots import make_subplots

import scvi

from cellbro.util.Param import *

scvi_model_params = ParamsDict([
    Param("gene_likelihood", "Gene Likelihood", default="zinb", type=list, allowed_values=["zinb", "nb", "poisson"]),
    Param("latent_distribution", "Latent Distribution", default="normal", type=list, allowed_values=["normal", "ln"]),
    Param("n_hidden", "Nodes Per Hidden Layer", default=128, type=int, _min=1, _max=1024, step=1),
    Param("n_latent", "Latent Dimension", default=10, type=int, _min=1, _max=1024, step=1),
    Param("n_layers", "Num. of Hidden Layers", default=1, type=int, _min=1, _max=10, step=1),
    Param("dropout_rate", "Dropout Rate", default=0.1, type=float, _min=0, _max=1, step=0.1),
])

scvi_train_params = ParamsDict([
    Param("max_epochs", "Max Epochs", default=1, type=int, _min=1, _max=10000, step=1),
    Param("train_size", "Train Size", default=0.9, type=float, _min=0, _max=1, step=0.1),
    Param("batch_size", "Batch Size", default=128, type=int, _min=1, _max=1024, step=1),
    # Param("early_stopping", "Early Stopping", default=True, type=bool),
])

def setup(dataset, params):
    scvi.model.SCVI.setup_anndata(dataset.adata, layer="counts", **params)

def fit(dataset, model_params, train_params):
    vae = scvi.model.SCVI(dataset.adata, **model_params)
    vae.train(**train_params)

    dataset.adata.obsm["X_scvi"] = vae.get_latent_representation()
    sc.pp.neighbors(dataset.adata, random_state=0, use_rep="X_scvi", key_added="neighbors_scvi")

    likelihood_params = vae.get_likelihood_parameters(dataset.adata)

    for key, value in likelihood_params.items():
        dataset.adata.layers[f"scvi_{key}"] = value
