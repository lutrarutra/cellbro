from enum import unique, Enum, auto

from flask import Flask, render_template, redirect, url_for
from sassutils.wsgi import SassMiddleware
from flask_sqlalchemy import SQLAlchemy

import scanpy as sc

adata = sc.read_h5ad("data/adata.h5ad")

db = SQLAlchemy()

class CellFeature(db.Model):
    key = db.Column(db.String(255), primary_key=True)
    feature_type = db.Column(db.SmallInteger(), nullable=False)

    def get_feature_type(self):
        return {
            0: "feature",
            1: "gene"
        }[self.feature_type]

from app import routes

def register_blueprints(app):
    app.register_blueprint(routes.api.projection_bp)
    app.register_blueprint(routes.api.query_bp)

def create_app():
    app = Flask(__name__)

    app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///:memory:"
    db.init_app(app)
    with app.app_context():
        db.create_all()

        for obs_key in adata.obs_keys():
            db.session.add(CellFeature(key=obs_key, feature_type=0))
        
        for var_name in sorted(adata.var_names):
            db.session.add(CellFeature(key=var_name, feature_type=1))
            
        db.session.commit()


    @app.route('/')
    def index_page():
        return redirect("/cells")

    @app.route('/cells')
    def cells_page():
        projection_types = [{"projection_type" : adata.obsm_keys()}]
        features = url_for("queries.querycellfeatures")

        return render_template(
            "cells.html", projection_types=projection_types, features=features,
            features_default=adata.obs_keys()[0]
        )
    
    @app.route('/genes')
    def genes_page():
        return render_template("genes.html")
    
    @app.route('/qc')
    def qc_page():
        return render_template("qc.html")
    
    @app.route('/files')
    def files_page():
        return render_template("files.html")
    
    register_blueprints(app)

    app.wsgi_app = SassMiddleware(app.wsgi_app, {
        "app" : ("static/sass", "static/css", "/static/css")
    })

    return app