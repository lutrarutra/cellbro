
import Window

class App():
    def __init__(self):
        self.window = Window.Window()

    def loop(self):
        while self.window.update():
            pass
