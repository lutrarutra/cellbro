import dash_bootstrap_components as dbc
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go

from cellbro.util.Param import *

import scout

figure_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
)


class PCA:
    def __init__(self, dataset, color, pc_x, pc_y, hist_n_pcs, hist_type):
        self.dataset = dataset
        self.color = color
        self.pc_x = pc_x - 1
        self.pc_y = pc_y - 1
        self.hist_n_pcs = hist_n_pcs
        self.hist_type = hist_type

        # if color in self.dataset.adata.obs_keys():
        #     self.color = self.dataset.adata.obs[color].values
        # else:
        #     self.color = (
        #         self.dataset.adata.X[:, self.dataset.adata.var.index.get_loc(color)]
        #         .toarray()
        #         .T[0]
        #     )
        # self.params = params

    def projection(self):
        fig = scout.ply.projection(
            self.dataset.adata, obsm_layer="X_pca", hue=self.color,
            layout=figure_layout, components=[self.pc_x, self.pc_y],
        )
        x_ratio = self.dataset.adata.uns["pca"]["variance_ratio"][self.pc_x]
        y_ratio = self.dataset.adata.uns["pca"]["variance_ratio"][self.pc_y]

        fig.update_layout(figure_layout)
        fig.update_xaxes(scaleanchor="y", scaleratio=1)
        fig.update_layout(
            dict(
                xaxis=dict(showticklabels=False),
                yaxis=dict(showticklabels=False),
            )
        )

        fig.update_layout(
            xaxis_title=f"PC {self.pc_x+1} ({x_ratio*100:.1f} %)",
            yaxis_title=f"PC {self.pc_y+1} ({y_ratio*100:.1f} %)",
        )

        return fig

    def explain_variance(self):
        if self.hist_type == "Bar":
            fig = px.bar(
                x=range(1, self.hist_n_pcs + 1),
                y=self.dataset.adata.uns["pca"]["variance_ratio"][: self.hist_n_pcs],
                labels={"x": f"PC", "y": f"Variance Ratio"},
            )
        else:
            y = self.dataset.adata.uns["pca"]["variance_ratio"][: self.hist_n_pcs]

            _plot = None
            if self.hist_type == "Line":
                _plot = px.line
            elif self.hist_type == "Area":
                _plot = px.area
            elif self.hist_type == "Cumulative":
                _plot = px.area
                y = np.cumsum(y)
            else:
                assert False, f"Unknown hist_type: {self.hist_type}"

            fig = _plot(
                x=range(1, self.hist_n_pcs + 1),
                y=y,
                labels={"x": f"PC", "y": f"Variance Ratio"},
                markers=True,
            )
            # fig.update_layout(yaxis_range=[-1, self.hist_n_pcs + 1])

        fig.update_layout(figure_layout)
        fig.update_layout(
            xaxis_title="PC",
            yaxis_title="Variance Ratio",
            margin=dict(t=10, b=10, l=10, r=10),
        )
        return fig

    def correlation_circle(self):
        x = self.dataset.adata.varm["PCs"][:, self.pc_x]
        y = self.dataset.adata.varm["PCs"][:, self.pc_y]
        x_ratio = self.dataset.adata.uns["pca"]["variance_ratio"][self.pc_x]
        y_ratio = self.dataset.adata.uns["pca"]["variance_ratio"][self.pc_y]

        dist = np.sqrt(x**2 + y**2)

        text = self.dataset.adata.var_names.values.copy()
        text[dist < np.quantile(dist, 0.999)] = ""

        df = pd.DataFrame(
            {"x": x, "y": y, "Gene": self.dataset.adata.var_names.values, "text": text}
        )
        df["textposition"] = ""
        df.loc[df["y"] > 0, "textposition"] = (
            df.loc[df["y"] > 0, "textposition"] + "top"
        )
        df.loc[df["y"] == 0, "textposition"] = (
            df.loc[df["y"] == 0, "textposition"] + "center"
        )
        df.loc[df["y"] < 0, "textposition"] = (
            df.loc[df["y"] < 0, "textposition"] + "bottom"
        )

        df.loc[df["x"] < 0, "textposition"] = (
            df.loc[df["x"] < 0, "textposition"] + " left"
        )
        df.loc[df["x"] > 0, "textposition"] = (
            df.loc[df["x"] > 0, "textposition"] + " right"
        )

        xmin, xmax = np.min(x), np.max(x)
        ymin, ymax = np.min(y), np.max(y)
        xaxis_range = [xmin - 0.3 * xmax, xmax + 0.3 * xmax]
        yaxis_range = [ymin - 0.3 * ymax, ymax + 0.3 * ymax]

        fig = px.scatter(
            df, x="x", y="y", text="text", hover_name="Gene",
            hover_data={"x": ":.2f", "y": ":.2f", "Gene": False, "text": False},
            labels={"x": f"PC {self.pc_x+1}", "y": f"PC {self.pc_y+1}"}
        )

        fig.update_traces(
            textposition=df["textposition"],
            marker=dict(size=5, line=dict(width=1, color="DarkSlateGrey"))
        )
        
        fig.update_layout(
            xaxis_title=f"PC {self.pc_x+1} ({x_ratio*100:.1f} %)",
            yaxis_title=f"PC {self.pc_y+1} ({y_ratio*100:.1f} %)",
            xaxis_range=xaxis_range, yaxis_range=yaxis_range,
        )
        fig.update_layout(figure_layout)
        return fig

    def explain_corr(self):
        cats = self.dataset.get_categoric()
        Rs = np.zeros((10, len(cats)))
        for i, cat in enumerate(cats):
            for j in range(10):
                Rs[j, i] = (
                    np.corrcoef(
                        self.dataset.adata.obs[cat].cat.codes,
                        self.dataset.adata.obsm["X_pca"][:, j],
                    )[0, 1]
                    ** 2
                )

        df = pd.DataFrame(Rs, columns=cats, index=range(1, 11)).reset_index()
        fig = px.scatter(df, x="index", y=cats)
        fig.update_traces(
            marker=dict(size=10, line=dict(width=2, color="DarkSlateGrey"))
        )
        fig.update_layout(
            xaxis=dict(
                title="PC", tick0=1, dtick=1, showgrid=False, zeroline=False,
                visible=True, showticklabels=True,
            ),
            yaxis=dict(
                range=[0, 1], title="R<sup>2</sup>", tick0=0, dtick=0.2, showgrid=False,
                zeroline=False, visible=True, showticklabels=True,
            ),
            legend_title="Feature",
        )
        return fig

    def bottom_plots(self):
        var_explained_fig = self.explain_variance()
        corr_explained_fig = self.explain_corr()

        return var_explained_fig, corr_explained_fig

    def plot(self):
        main_plot = self.projection()
        secondary_plot = self.correlation_circle()
        var_explained_fig, corr_explained_fig = self.bottom_plots()
        return [main_plot, secondary_plot, var_explained_fig, corr_explained_fig]

        # for param in pca_params.values():
        #     inputs[param.key] = Input(component_id=f"pca-{param.key}", component_property="value")
        return inputs
