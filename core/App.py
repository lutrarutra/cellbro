
import Window, util.Logger as Logger
import util.Dataset as Dataset
from analysis import preprocessing as pp

class App():
    def __init__(self):
        self.logger = Logger.Logger()
        self.window = Window.Window(self, self.logger)
        self.dataset = None
        self.processes = {}

    def loop(self):
        while self.window.update():
            pass

    def init_dataset(self, path, file_type):
        self.dataset = Dataset.Dataset(path, file_type, self.logger)
        self.window.childs["preprocessing_progress"] = pp.PreprocessProgress()
        self.window.childs["ask_phenodata"] = pp.AnnotateForm()
        
