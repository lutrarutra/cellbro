from abc import ABC, abstractmethod

import dash_bootstrap_components as dbc
from dash import Input, Output, State, dcc, html, ctx

from cellbro.util.Param import ParamsDict, Param

class ProjectionType(ABC):
    def __init__(self, dataset, params: ParamsDict):
        self.dataset = dataset
        self.params = params

    @classmethod
    @abstractmethod
    def get_key(cls) -> str:
        ...

    @staticmethod
    @abstractmethod
    def get_params() -> ParamsDict:
        ...

    @abstractmethod
    def apply(self) -> str:
        ...

    @classmethod
    def parse_params(cls, params):
        key = cls.get_key()
        projection_params = dict(
            [
                (param_key.replace(f"{key}_", ""), params[param_key])
                for param_key in params.keys()
                if param_key.startswith(f"{key}_")
            ]
        )
        return projection_params

    # def add_params(self):
    #     if f"{self.get_type()}_params" not in self.dataset.adata.uns.keys():
    #         self.dataset.adata.uns[f"{self.get_type()}_params"] = {}

    #     rerun = (
    #         self.dataset.adata.uns[f"{self.get_type()}_params"] != self.params.unravel()
    #     )
    #     self.dataset.adata.uns[f"{self.get_type()}_params"] = self.params.unravel()

    #     return rerun

    @classmethod
    def get_layout(cls, page_id_prefix):
        divs = []
        for key, param in cls.get_params().items():
            if param.type == bool:
                inp = dbc.Switch(
                    id=f"{page_id_prefix}-projection-{cls.get_key()}-{key}",
                    value=param.default,
                )
            else:
                inp = dbc.Input(
                    id=f"{page_id_prefix}-projection-{cls.get_key()}-{key}",
                    type=param.input_type,
                    value=param.value,
                    step=param.step if param.step != None else 0.1,
                    className="param-input",
                    placeholder=param.placeholder,
                )

            divs.append(
                html.Div(
                    children=[
                        html.Label(
                            param.name,
                            className="param-label",
                        ),
                        inp,
                    ],
                    className="param-row",
                )
            )

        return divs
