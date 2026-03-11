import warnings
from .MessageRedis import SyncMessageRedis

class WarningCaptureHandler:
    def __init__(self, task_id: str, message_client: SyncMessageRedis):
        self._task_id = task_id
        self.client = message_client
        self._channel = "stdout"
        self._original_showwarning = None
        
    def custom_showwarning(self, message, category, filename, lineno, file=None, line=None):
        warning_msg = warnings.formatwarning(message, category, filename, lineno, line)
        self.client.warning(self._task_id,  warning_msg)

    def start(self):
        self._original_showwarning = warnings.showwarning
        warnings.showwarning = self.custom_showwarning
    
    def stop(self):
        if self._original_showwarning:
            warnings.showwarning = self._original_showwarning
    
    def __enter__(self):
        self.start()
        return self
    
    def __exit__(self, exc_type, exc_value, traceback):
        self.stop()