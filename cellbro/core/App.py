import tkinter as tk
from tkinter import filedialog

from dash import Dash, html, dcc
import plotly.express as px
import plotly.graph_objects as go

from cellbro.util.Dataset import Dataset
from cellbro.core.Dashboard import Dashboard

class App():
    def __init__(self):
        self.dash_app = Dash(__name__)
        root = tk.Tk()
        root.withdraw()
        path = filedialog.askopenfilename()
        self.dataset = Dataset(path)
        self.dashboard = Dashboard(self)

    def run(self):
        self.dash_app.layout = self.dashboard.create_layout()
        self.dash_app.run_server(debug=False, host="127.0.0.1")