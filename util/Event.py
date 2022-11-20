import imgui

class Event():
    def __init__(self, key, callback=None, args={}):
        self.key = key
        self.args = args
        self.callback = callback

class EventHandler():
    def __init__(self, app):
        self.app = app
        self.events = {}
        self.completed_events = []

    def add_event(self, key, callback=None, args={}):
        self.events[key] = Event(key, callback, args=args)

    def complete_event(self, event_key, args=None):
        print("Event completed: " + event_key)
        if args is not None:
            self.events[event_key].args.update(args)
        self.completed_events.append(event_key)

    def remove_event(self, event_key):
        print("Event removed: " + event_key)
        del self.events[event_key]
