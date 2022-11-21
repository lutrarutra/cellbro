import random, string, multiprocessing

import scanpy as sc

class Figure():
    def __init__(self, app, type, plot_fnc):
        self.app = app
        self.id = f"{type}_" + "".join(random.choice(string.ascii_uppercase + string.digits) for _ in range(5))
        self.plot_fnc = plot_fnc

    def plot(self, kwargs):
        self.app.figures[self.id] = multiprocessing.Process(target=self.plot_fnc, kwargs=kwargs)
        self.app.figures[self.id].start()

