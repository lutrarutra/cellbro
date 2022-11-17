from time import gmtime, strftime
from enum import Enum

class LogLevel(Enum):
    DEBUG = 0
    INFO = 1
    WARNING = 2
    ERROR = 3

class Logger():
    def __init__(self, log_level=LogLevel.DEBUG):
        self.log = []
        self.log_level = log_level

    def _log(self, level, message):
        self.log.append((level, message, strftime("%H:%M", gmtime())))

    def debug(self, message):
        self._log(LogLevel.DEBUG, message)

    def info(self, message):
        self._log(LogLevel.INFO, message)

    def warning(self, message):
        self._log(LogLevel.WARNING, message)

    def error(self, message):
        self._log(LogLevel.ERROR, message)