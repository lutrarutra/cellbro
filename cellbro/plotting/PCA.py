import plotly.graph_objects as go
import plotly.express as px
from dash import html, dcc, Input, Output, State
import dash_bootstrap_components as dbc
from plotly.subplots import make_subplots

from cellbro.util.Param import *

import scanpy as sc
import pandas as pd
import numpy as np

figure_layout = go.Layout(
    paper_bgcolor='white',
    plot_bgcolor='white',
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=0, b=0, l=0, r=0),
)

pca_params = ParamsDict([
#     Param(key="pc_1", name="PC 1", default=1, type=int, description="", step=1, _min=0),
#     Param(key="pc_2", name="PC 2", default=2, type=int, description="", step=1, _min=0),
])

class PCA():
    def __init__(self, dataset, color, pc_x, pc_y, hist_n_pcs, hist_type):
        self.dataset = dataset
        self.color = color
        self.pc_x = pc_x-1
        self.pc_y = pc_y-1
        self.color_label = color
        self.hist_n_pcs = hist_n_pcs
        self.hist_type = hist_type

        if color in self.dataset.adata.obs_keys():
            self.color = self.dataset.adata.obs[color].values
        else:
            self.color = self.dataset.adata.X[:,self.dataset.adata.var.index.get_loc(color)].toarray().T[0]
        # self.params = params

    def projection(self):
        x_ratio = self.dataset.adata.uns["pca"]["variance_ratio"][self.pc_x]
        y_ratio = self.dataset.adata.uns["pca"]["variance_ratio"][self.pc_y]
        fig = px.scatter(
            x=self.dataset.adata.obsm["X_pca"][:,self.pc_x],
            y=self.dataset.adata.obsm["X_pca"][:,self.pc_y],
            color=self.color, color_discrete_sequence=sc.pl.palettes.default_20,
            labels={"x": f"PC {self.pc_x+1}", "y": f"PC {self.pc_y+1}"}
        )
        figure_layout["legend"] = dict(
            title=self.color_label.capitalize()
        )
        fig.update_layout(figure_layout)
        fig.update_xaxes(scaleanchor="y", scaleratio=1)
        fig.update_layout(dict(
            xaxis=dict(showticklabels=False),
            yaxis=dict(showticklabels=False),
        ))

        fig.update_layout(xaxis_title=f"PC {self.pc_x+1} ({x_ratio*100:.1f} %)", yaxis_title=f"PC {self.pc_y+1} ({y_ratio*100:.1f} %)")

        return fig

    def explain_variance(self):
        if self.hist_type == "Bar":
            fig = px.bar(
                x=range(1, self.hist_n_pcs + 1),
                y=self.dataset.adata.uns["pca"]["variance_ratio"][:self.hist_n_pcs],
                labels={"x": f"PC", "y": f"Variance Ratio"}
            )
        else:
            y = self.dataset.adata.uns["pca"]["variance_ratio"][:self.hist_n_pcs]

            if self.hist_type == "Line":
                _plot = px.line
            elif self.hist_type == "Area":
                _plot = px.area
            elif self.hist_type == "Cumulative":
                _plot = px.area
                y = np.cumsum(y)

            fig = _plot(
                x=range(1, self.hist_n_pcs + 1), y=y,
                labels={"x": f"PC", "y": f"Variance Ratio"}, markers=True
            )
            # fig.update_layout(yaxis_range=[-1, self.hist_n_pcs + 1])
        
        fig.update_layout(figure_layout)
        fig.update_layout(
            xaxis_title="PC", yaxis_title="Variance Ratio",
            margin=dict(t=10, b=10, l=10, r=10)
        )
        return fig

    def correlation_circle(self):
        x = self.dataset.adata.varm["PCs"][:,self.pc_x]
        y = self.dataset.adata.varm["PCs"][:,self.pc_y]
        x_ratio = self.dataset.adata.uns["pca"]["variance_ratio"][self.pc_x]
        y_ratio = self.dataset.adata.uns["pca"]["variance_ratio"][self.pc_y]

        dist = np.sqrt(x**2 + y**2)

        text = self.dataset.adata.var_names.values.copy()
        text[dist < np.quantile(dist, 0.999)] = ""



        df = pd.DataFrame({"x":x, "y":y, "Gene":self.dataset.adata.var_names.values,"text":text})
        df["textposition"] = ""
        df.loc[df["y"] > 0, "textposition"] = df.loc[df["y"] > 0, "textposition"] + "top"
        df.loc[df["y"] == 0, "textposition"] = df.loc[df["y"] == 0, "textposition"] + "center"
        df.loc[df["y"] < 0, "textposition"] = df.loc[df["y"] < 0, "textposition"] + "bottom"

        df.loc[df["x"] < 0, "textposition"] = df.loc[df["x"] < 0, "textposition"] + " left"
        df.loc[df["x"] > 0, "textposition"] = df.loc[df["x"] > 0, "textposition"] + " right"

        # df.loc[(df["x"] > 0) & (df["y"] > 0), "textposition"] = "top right"
        # df.loc[(df["x"] > 0) & (df["y"] < 0), "textposition"] = "top left"
        # df.loc[(df["x"] < 0) & (df["y"] > 0), "textposition"] = "bottom right"
        # df.loc[(df["x"] < 0) & (df["y"] < 0), "textposition"] = "left bottom"
        
        xmin, xmax = np.min(x), np.max(x)
        ymin, ymax = np.min(y), np.max(y)
        xaxis_range = [xmin - 0.3 * xmax, xmax + 0.3 * xmax]
        yaxis_range = [ymin - 0.3 * ymax, ymax + 0.3 * ymax]

        fig = px.scatter(
            df, x="x", y="y", text="text", hover_name="Gene",
            hover_data={"x":":.2f", "y":":.2f", "Gene":False, "text":False},
            labels={"x": f"PC {self.pc_x+1}", "y": f"PC {self.pc_y+1}"}
            # color_discrete_sequence=sc.pl.palettes.default_20,
        )

        fig.update_traces(textposition=df["textposition"])
        fig.update_layout(
            xaxis_title=f"PC {self.pc_x+1} ({x_ratio*100:.1f} %)",
            yaxis_title=f"PC {self.pc_y+1} ({y_ratio*100:.1f} %)",
            xaxis_range=xaxis_range, yaxis_range=yaxis_range,
        )
        fig.update_layout(figure_layout)
        return fig

    def bottom_plots(self):
        fig = make_subplots(rows=1, cols=3)
        fig.add_trace(
            list(self.explain_variance().select_traces())[0],
        )
        fig.update_layout(figure_layout)
        return fig

    def plot(self):
        main_plot = self.projection()
        secondary_plot = self.correlation_circle()
        bottom_plots = self.bottom_plots()
        return [main_plot, secondary_plot, bottom_plots]
    
    @staticmethod
    def get_callback_outputs():
        return [
            Output(component_id="pca-main-plot", component_property="figure"),
            Output(component_id="pca-secondary-plot", component_property="figure"),
            Output(component_id="pca-bottom-plot", component_property="figure"),
        ]
        
    # Inputs to Projection
    @staticmethod
    def get_callback_inputs():
        inputs = {
            "color": Input(component_id="pca-color", component_property="value"),
            "pc_x": Input(component_id="pca-x-component", component_property="value"),
            "pc_y": Input(component_id="pca-y-component", component_property="value"),
            "hist_type": Input(component_id="pca-hist_type", component_property="value"),
            "hist_n_pcs": Input(component_id="pca-hist_n_pcs", component_property="value"),
        }

        # for param in pca_params.values():
        #     inputs[param.key] = Input(component_id=f"pca-{param.key}", component_property="value")
        return inputs

    @staticmethod
    def create_layout(dataset):
        top_sidebar = html.Div(children=[
            html.Div([
                html.H3("PCA Projection Settings"),
            ], className="sidebar-header"),
            dcc.Loading(type="circle", children=[
                html.Div(children=[
                    PCA.params_layout(),
                ], className="sidebar-parameters"),
                html.Div([
                    dbc.Button("Filter", color="primary", className="mr-1", id="pca-submit"),
                ], className="sidebar-footer")
            ],),
        ], className="top-sidebar sidebar")


        bottom_sidebar = html.Div(children=[
            html.Div([
                html.H3("PCA Plots"),
            ], className="sidebar-header"),
            dcc.Loading(type="circle", children=[
                html.Div(children=[
                    html.Label("Plot Type"),
                    dcc.Dropdown(["Bar", "Line", "Area", "Cumulative"], value="Bar", id="pca-hist_type", clearable=False),
                ], className="param-row"),
                # Num components
                html.Div(children=[
                    html.Label("Num. Components"),
                    dcc.Input(
                        id="pca-hist_n_pcs", type="number", value=30, min=2, step=1,
                        max=dataset.adata.uns["pca"]["variance_ratio"].shape[0]+1,
                        className="param-input"
                    ),
                ], className="param-row"),
            ],)
        ], id="pca-bottom-sidebar", className="bottom-sidebar sidebar")

        main_figure = html.Div(children=[
            html.Div(children=[
                # Projection Color
                html.Div(children=[
                    html.Label("Color"),
                    dcc.Dropdown(dataset.adata.obs_keys() + dataset.adata.var_names.tolist(), value=dataset.adata.obs_keys()[0], id="pca-color", clearable=False),
                ], className="param-column"),

                # X-axis component
                html.Div(children=[
                    html.Label("X Component"),
                    dcc.Input(
                        id="pca-x-component", type="number", value=1, min=1, step=1,
                        max=dataset.adata.uns["pca"]["variance_ratio"].shape[0]+1,
                        className="param-input"
                    ),
                ], className="param-column"),

                # X-axis component
                html.Div(children=[
                    html.Label("Y Component"),
                    dcc.Input(
                        id="pca-y-component", type="number", value=2, min=1, step=1,
                        max=dataset.adata.uns["pca"]["variance_ratio"].shape[0]+1,
                        className="param-input"
                    ),
                ], className="param-column"),

            ], id="pca-main-select", className="top-parameters"),

            html.Div([
                dcc.Loading(
                    id="loading-pca-main", type="circle",
                    children=[html.Div(dcc.Graph(id="pca-main-plot", className="main-plot"))],
                )
            ], id="pca-main-figure", className="main-figure")
        ], className="main")

        secondary_figure = html.Div(children=[
            html.Div(children=[
            ], id="pca-secondary-select", className="top-parameters"),
            html.Div([
                dcc.Loading(
                    id="loading-pca-secondary", type="circle",
                    children=[html.Div(dcc.Graph(id="pca-secondary-plot", className="secondary-plot"))],
                )
            ], id="pca-secondary-figure", className="secondary-figure")
        ], className="secondary")

        bottom_figure = html.Div(children=[
            dcc.Loading(
                id="loading-pca-bottom", className="loading-bottom", type="circle",
                children=[html.Div(dcc.Graph(id="pca-bottom-plot", className="bottom-plot"))],
            )
        ], id="pca-bottom-figure", className="bottom-figure")

        return top_sidebar, main_figure, secondary_figure, bottom_sidebar, bottom_figure

    @staticmethod
    def params_layout():
        divs = []
        for key, param in pca_params.items():
            divs.append(html.Div(children=[
                html.Label(
                    param.name, className="param-label",
                ),
                dcc.Input(
                    id=f"pca-{key}", type=param.input_type, value=param.value,
                    step=param.step if param.step != None else 0.1, className="param-input",
                ),
            ], className="param-row"))

        layout = html.Div(children=[divs])
        return layout