from functools import cmp_to_key

import plotly.graph_objects as go

import scout

from ...util.Param import Param, ParamsDict

default_layout = go.Layout(
    paper_bgcolor="white",
    plot_bgcolor="white",
    xaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    yaxis=dict(showgrid=False, zeroline=False, visible=True, showticklabels=True),
    margin=dict(t=5, b=5, l=5, r=5),
)

de_params = ParamsDict(
    [
        Param(
            "method",
            "Method",
            default="t-test",
            type=str,
            description="",
            allowed_values=["logreg", "t-test", "wilcoxon", "t-test_overestim_var"],
        ),
        Param(
            "corr_method",
            "Correction Method",
            default="benjamini-hochberg",
            type=str,
            description="",
            allowed_values=[
                "benjamini-hochberg",
                "bonferroni",
            ],
        ),
    ]
)


# def apply(dataset, params):
#     if params["submit"] is None:
#         groupby_options = get_groupby_options(dataset=dataset)
#         if len(groupby_options) == 0:
#             return dict(update=False)

#         groupby_selected = groupby_options[0]
#         refs_options = get_reference_options(dataset=dataset, groupby=groupby_selected)
#         refs_selected = next(
#             iter(dataset.adata.uns[f"rank_genes_{groupby_selected}"].keys())
#         )
        
#         return dict(update=True)

#     key = f"rank_genes_{params['groupby']}"

#     if key not in dataset.adata.uns.keys():
#         scout.tl.rank_marker_genes(dataset.adata, groupby=params["groupby"])

#     groupby_options = get_groupby_options(dataset=dataset)
#     refs_options = get_reference_options(dataset, params["groupby"])

#     groupby_selected = params["groupby"]
#     refs_selected = next(iter(dataset.adata.uns[key].keys()))

#     return dict(update=True)


def get_reference_options(dataset, groupby, target="Rest"):
    groupby = dataset.adata.uns.get(f"rank_genes_{groupby}", None)

    def _compare(a, b):
        if a < b:
            return -1
        elif a > b:
            return 1

        return 0

    def compare(item1, item2):
        x_words = item1.split(" vs. ")
        y_words = item2.split(" vs. ")
        try:
            xt = int(x_words[0])
            yt = int(y_words[0])
            return _compare(xt, yt)
        except:
            if x_words[0] == y_words[0]:
                try:
                    xt = int(x_words[1])
                    yt = int(y_words[1])
                    return _compare(xt, yt)
                except:
                    return _compare(x_words[1], y_words[1])
            else:
                return _compare(x_words[0], y_words[0])
        

    return sorted(list(groupby.keys()), key=cmp_to_key(compare))
