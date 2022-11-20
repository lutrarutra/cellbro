import imgui

from datetime import datetime
from enum import Enum
from util import gui

class LogLevel(Enum):
    DEBUG = 0
    INFO = 1
    WARNING = 2
    ERROR = 3

class Logger():
    def __init__(self, log_level=LogLevel.DEBUG):
        self.log = []
        self.log_level = log_level
        self.flags = gui.compute_flags([
            imgui.WINDOW_NO_RESIZE,
            imgui.WINDOW_NO_MOVE,
        ])

    def draw(self):
        if imgui.is_window_collapsed():
            return
        
        for level, message, time in self.log:
            if level.value >= self.log_level.value:
                if level.value == LogLevel.DEBUG.value:
                    imgui.push_style_color(imgui.COLOR_TEXT, 0, 0.7, 0.7)
                elif level.value == LogLevel.INFO.value:
                    imgui.push_style_color(imgui.COLOR_TEXT, 0.5, 0.5, 0.5)
                elif level.value == LogLevel.WARNING.value:
                    imgui.push_style_color(imgui.COLOR_TEXT, 1.0, 0.7, 0.0)
                elif level.value == LogLevel.ERROR.value:
                    imgui.push_style_color(imgui.COLOR_TEXT, 1.0, 0.0, 0.0)

                imgui.text(f"[{time}] {message}")
                imgui.pop_style_color()


    def _log(self, level, message):
        imgui.set_scroll_y(imgui.get_scroll_max_y())
        self.log.append((level, message, datetime.now().strftime("%H:%M")))

    def debug(self, message):
        self._log(LogLevel.DEBUG, message)

    def info(self, message):
        self._log(LogLevel.INFO, message)

    def warning(self, message):
        self._log(LogLevel.WARNING, message)

    def error(self, message):
        self._log(LogLevel.ERROR, message)