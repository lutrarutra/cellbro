class Form():
    def __init__(self, event_handler, event_key):
        self.finished = False
        self.event_handler = event_handler
        self.event_key = event_key