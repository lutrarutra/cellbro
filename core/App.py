
import Window, util.Logger as Logger
import util.Dataset as Dataset

class App():
    def __init__(self):
        self.logger = Logger.Logger()
        self.window = Window.Window(self, self.logger)
        self.dataset = None

    def loop(self):
        while self.window.update():
            pass

    def init_dataset(self, path, file_type):
        self.dataset = Dataset.Dataset(path, file_type, self.logger)
