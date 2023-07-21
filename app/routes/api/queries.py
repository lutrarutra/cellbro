import json
from flask_restful import Resource, Api
from flask import Blueprint, request

query_bp = Blueprint("queries", __name__)
query_api = Api(query_bp)

from app import CellFeature

class QueryCellFeatures(Resource):
    def get(self):
        word = request.args.get("query")

        results = CellFeature.query.filter(
            CellFeature.key.like("%{}%".format(word))
        ).order_by(CellFeature.feature_type).limit(10).all()

        for result in results:
            result.feature_type = result.get_feature_type()

        return [{"key": result.key, "feature_type": result.feature_type} for result in results]
    
query_api.add_resource(QueryCellFeatures, "/cell_features")
