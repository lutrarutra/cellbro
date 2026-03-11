import sys
from io import StringIO
import logging

from .MessageRedis import SyncMessageRedis


class StdoutCaptureHandler(logging.Handler):
    def __init__(self, task_id: str, client: SyncMessageRedis):
        super().__init__()
        self._original_stdout = sys.stdout
        self._stdout = StringIO()
        self._task_id = task_id
        self.client =client
        self._channel = "stdout"
    
    def write(self, text):
        """Handle print() calls when sys.stdout = self"""
        if not text or text == '\n':
            return
        text = text.rstrip('\n')
        self._stdout.write(text + '\n')
        self.client.log(self._task_id, text, category="log")
        
    def emit(self, record):
        msg = self.format(record)
        self._stdout.write(msg + '\n')
        self.client.log(self._task_id, msg, category="log")
    
    def start(self):
        sys.stdout = self
    
    def stop(self):
        sys.stdout = self._original_stdout
    
    def get_captured_output(self):
        return self._stdout.getvalue()
    
    def __enter__(self):
        self.start()
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.stop()