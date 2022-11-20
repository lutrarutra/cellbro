import threading

class TaskHandler():
    def __init__(self):
        self.tasks = {}

    def add_task(self, key, task):
        self.tasks[key] = task

    def is_blocking(self):
        for task_key in self.tasks:
            if self.tasks[task_key].blocking and self.tasks[task_key].is_running():
                return True
        return False

    def process_tasks(self):
        for task_key in list(self.tasks.keys()):
            if not self.tasks[task_key].is_running():
                print("Task completed: " + task_key)
                self.tasks[task_key].join()
                del self.tasks[task_key]
    

class Task():
    def __init__(self, target, args=(), kwargs={}, blocking=True, on_finish=None):
        self.target = target
        self.args = args
        self.blocking = blocking
        self.process = threading.Thread(target=self.target, args=self.args, kwargs=kwargs)
        self.on_finish = on_finish

    def is_running(self):
        return self.process.is_alive()

    def start(self):
        self.process.start()

    def complete(self):
        print("Task completed: " + self.key)
        self.finished = True

    def join(self):
        self.process.join()
