import json
from flask_restful import Resource, Api
from flask import Blueprint, request, session

projection_bp = Blueprint("projection", __name__, url_prefix="/api/plotting")
projection_api = Api(projection_bp)

import plotly
from cellbro.plotting import projection

from cellbro.app import adata, logger

class Projection(Resource):
    def post(self):
        content = request.get_json()

        fig = projection.projection(
            adata, color=content["color"], obsm_layer=content["type"],
            width=content["width"], height=content["height"] 
        )

        session["selected_projection_type"] = content["type"]
        session["selected_feature"] = content["color"]

        logger.debug(f"Projection plot: color='{content['color']}', type='{content['type']}'")

        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    
    
projection_api.add_resource(Projection, "/projection")
