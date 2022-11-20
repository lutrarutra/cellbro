import imgui

import scanpy as sc

from util import Task

class Leiden():
    def __init__(self, app):
        self.app = app
        self.params = {
            "resolution": 1.0,
            "random_state": 0,
        }

    def draw(self):
        for key, value in self.params.items():
            if isinstance(value, int):
                _, self.params[key] = imgui.input_int(key, value)
            elif isinstance(value, float):
                _, self.params[key] = imgui.input_float(key, value)
        
        imgui.set_cursor_pos((imgui.get_cursor_pos()[0], imgui.get_window_height() - 40))

    def apply(self, start_thread=False):
        if start_thread:
            self.app.task_handler.add_task(
                "leiden_clustering",
                Task.Task(target=sc.tl.leiden, args=(self.app.dataset.adata,), kwargs=self.params)
            )
            self.app.task_handler.tasks["leiden_clustering"].start()
        else:
            sc.tl.leiden(self.app.dataset.adata, **self.params)