# import dash_bootstrap_components as dbc
# import numpy as np
# import pandas as pd
# import plotly.express as px
import plotly.graph_objects as go

# from cellbro.util.Param import *

# import scout

default_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
)