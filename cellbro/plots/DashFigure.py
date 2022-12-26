from abc import ABC, abstractmethod
from typing import Literal

from cellbro.util.Components import DashComponent

class DashFigure(DashComponent, ABC):
    def __init__(self, dataset, page_id_prefix, loc_class: Literal["main", "secondary", "bottom"]):
        super().__init__(page_id_prefix)
        self.dataset = dataset
        self.loc_class = loc_class

    @abstractmethod
    def create_layout(self) -> list:
        ...

    def get_sidebar_params(self) -> list:
        return []

