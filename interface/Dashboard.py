import imgui

class Dashboard():
    def __init__(self, window):
        self.window = window
        self.parts = {
            "log":None
        }


    def draw(self):
        for part in self.parts:
            self.parts[part].draw()