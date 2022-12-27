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
    def __init__(self, dataset, page_id_prefix):
        self.dataset = dataset
        self.page_id_prefix = page_id_prefix

    # @abstractmethod
    # def create_layout(self) -> list:
    #     ...

    @abstractmethod
    def setup_callbacks(self, app):
        ...

class PlotAction(DashAction, ABC):
    def __init__(self, dataset, page_id_prefix, loc_class):
        super().__init__(dataset, page_id_prefix)
        self.loc_class = loc_class
