import sys
from io import StringIO
import logging
from celery.utils.log import get_task_logger

logger = get_task_logger(__name__)

class StdoutCaptureHandler(logging.Handler):
    def __init__(self, task_id: str, redis_client):
        super().__init__()
        self._original_stdout = sys.stdout
        self._stdout = StringIO()
        self._task_id = task_id
        self.redis_client = redis_client
        self._channel = f"task:{task_id}"
    
    def write(self, text):
        """Handle print() calls when sys.stdout = self"""
        if not text or text == '\n':
            return
        text = text.rstrip('\n')
        self._stdout.write(text + '\n')
        self.redis_client.publish(self._channel, text)
        
    def emit(self, record):
        msg = self.format(record)
        self._stdout.write(msg + '\n')
        self.redis_client.publish(self._channel, msg)
    
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