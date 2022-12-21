from abc import ABC, abstractmethod

import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html
from dash.exceptions import PreventUpdate
from plotly.subplots import make_subplots

from cellbro.util.Dataset import Dataset
from cellbro.util.Param import *
# from cellbro.core.app import App


class DashAction(ABC):
    # params: ParamsDict
    # callbacks: list
    # dataset: Dataset

    def __init__(self, dataset):
        self.dataset = dataset

    # @abstractmethod
    # def create_layout(self) -> list:
    #     ...

    @abstractmethod
    def setup_callbacks(self, app):
        ...

    def get_dummies(self):
        return []
