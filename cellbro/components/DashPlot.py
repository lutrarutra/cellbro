from abc import ABC, abstractmethod
from typing import Literal

from .DashComponent import DashComponent

class DashPlot(DashComponent, ABC):
    def __init__(self, dataset, page_id, loc_class):
        super().__init__(page_id, loc_class, "plot")
        self.dataset = dataset

    @abstractmethod
    def create_layout(self):
        ...

    def get_sidebar_params(self) -> list:
        return []

