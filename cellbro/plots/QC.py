import dash_bootstrap_components as dbc
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
import scanpy as sc
from dash import Input, Output, State, dcc, html
from plotly.subplots import make_subplots

import cellbro.util.Components as Components
from cellbro.util.Param import *

figure_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=0, b=0, l=0, r=0),
)

qc_params = ParamsDict(
    [
        Param(
            key="pct_counts_mt",
            name="Max MT%",
            default=5.0,
            type=float,
            description="",
            step=1,
        ),
        Param(
            key="min_genes",
            name="Min. Genes (per cell)",
            default=200,
            type=int,
            description="",
            step=100,
        ),
        Param(
            key="min_cells",
            name="Min. Cells (per gene)",
            default=10,
            type=int,
            description="",
            step=10,
        ),
    ]
)

def filter(self, submit):
    # Makes sure that filtering is not done on initial load
    if submit is None:
        return list(self.params.values()) # !=

    self.dataset.adata = self.dataset.adata[
        self.dataset.adata.obs.pct_counts_mt < self.params["pct_counts_mt"], :
    ].copy()
    sc.pp.filter_cells(self.dataset.adata, min_genes=self.params["min_genes"])
    sc.pp.filter_genes(self.dataset.adata, min_cells=self.params["min_cells"])

    return list(self.params.values())

def dispersion_plot(dataset, params):
    fig = px.scatter(
        dataset.adata.var.reset_index(),
        x="mu",
        y="cv2",
        log_x=True,
        log_y=True,
        color_continuous_scale=px.colors.sequential.Viridis,
        hover_name="index",
    )
    fig.update_traces(
        marker=dict(size=5, line=dict(width=1, color="DarkSlateGrey"))
    )
    fig.update_layout(figure_layout)
    fig.update_layout(xaxis_title="Log Mean Expression", yaxis_title="CV^2")
    return fig

def violin_plot(dataset, params):
    violins = ["n_genes_by_counts", "total_counts", "pct_counts_mt"]
    fig = make_subplots(rows=1, cols=3)
    x = dataset.adata.obs[violins]
    df = pd.DataFrame(x, columns=violins)

    for i, feature in enumerate(violins):
        fig.add_trace(
            go.Violin(
                name=feature.replace("_", " ").title(),
                y=df[feature],
                box_visible=True,
                points="all",
                pointpos=0,
                marker=dict(size=2),
                jitter=0.6,
            ),
            row=1,
            col=i + 1,
        )

    fig.update_layout(figure_layout)
    fig.update_layout(showlegend=False)

    return fig

def mt_plot(dataset, params):
    color = (
        dataset.adata.obs["pct_counts_mt"] > params["pct_counts_mt"]
    ).values

    cmap = ["#d3d3d3", "#636EFA"] if color[0] else ["#636EFA", "#d3d3d3"]

    fig = px.scatter(
        dataset.adata.obs,
        x="total_counts",
        y="pct_counts_mt",
        color=color,
        color_discrete_sequence=cmap,
    )
    fig.update_traces(
        marker=dict(size=5, line=dict(width=1, color="DarkSlateGrey"))
    )
    fig.update_layout(figure_layout)
    fig.update_layout(
        xaxis_title="Total Counts",
        yaxis_title="MT%",
        legend_title_text=f"MT > {params['pct_counts_mt']:.1f}%",
    )
    fig.add_hline(
        y=params["pct_counts_mt"],
        line_width=1,
        line_dash="dash",
        line_color=sc.pl.palettes.default_20[3],
    )
    return fig

def plot(dataset, params):
    main_figure = mt_plot(dataset, params)
    bottom_figure = violin_plot(dataset, params)
    secondary_figure = dispersion_plot(dataset, params)

    return [main_figure, secondary_figure, bottom_figure]


def params_layout():
    layout = html.Div(children=[
        html.Div(children=[
            html.Label(
                param.name, className="param-label",
            ),
            dcc.Input(
                id=f"qc-{key}", type=param.input_type, value=param.value,
                step=param.step if param.step != None else 0.1, className="param-input",
            ),
        ], className="param-row") for key, param in qc_params.items()
    ])
    return layout