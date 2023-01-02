from dash import dcc

from .DashComponent import DashComponent

class DashStore(DashComponent):
    def __init__(self, page_id, loc_class, type, storage_type="local"):
        super().__init__(page_id, loc_class, type)
        self.storage_type = storage_type

    def create_layout(self):
        return dcc.Store(id=self.cid.to_dict(), storage_type=self.storage_type)
