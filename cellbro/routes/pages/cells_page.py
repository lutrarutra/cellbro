import json
from flask_restful import Resource, Api
from flask import Blueprint, request, session, render_template, url_for

from cellbro.app import adata, logger

def render_cells_page():
    projection_types = [{"projection_type" : adata.obsm_keys()}]
    features = url_for("queries.querycellfeatures")

    if "selected_projection_type" in session.keys():
        selected_projection_type = session["selected_projection_type"]
    else:
        selected_projection_type = projection_types[0]["projection_type"][0]

    if "selected_feature" in session.keys():
        selected_feature = session["selected_feature"]
    else:
        selected_feature = adata.obs_keys()[0]

    return render_template(
        "cells.html", projection_types=projection_types, features=features,
        selected_projection_type=selected_projection_type,
        selected_feature=selected_feature
    )