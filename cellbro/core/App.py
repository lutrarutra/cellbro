# Run this app with `python app.py` and
# visit http://127.0.0.1:8050/ in your web browser.

from dash import Dash, html, dcc
import plotly.express as px
import plotly.graph_objects as go

class App():
    def __init__(self):
        self.app = Dash(__name__)
        self.app.run_server(debug=True, host="127.0.0.1")
    