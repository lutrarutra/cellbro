import json
from flask_restful import Resource, Api
from flask import Blueprint, request

projection_bp = Blueprint("projection", __name__)
projection_api = Api(projection_bp)

import plotly
from app.plotting import projection

from app import adata

class Projection(Resource):
    def post(self):
        content = request.get_json()

        fig = projection.projection(
            adata, color=content["color"], obsm_layer=content["type"],
            width=content["width"], height=content["height"] 
        )

        print(content["color"])

        return json.dumps(fig, cls=plotly.utils.PlotlyJSONEncoder)
    
    
projection_api.add_resource(Projection, "/projection")
