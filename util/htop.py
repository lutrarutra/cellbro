import psutil, time, array

import imgui

psutil.virtual_memory()

class HTOP():
    def __init__(self, frequency=1, history=100):
        self.n_threads = psutil.cpu_count()

        self.cpu_util = [0.0] * self.n_threads
        self.cpu_util_history = [0.0] * history

        self.memory_total_mb = psutil.virtual_memory().total / (1024 * 1024)
        self.memory_util = 0.0
        self.memory_available_mb = 0.0
        self.memory_util_history = [0.0] * history

        self.fps_history = [0.0]*10
        self.fps = 0.0

        self.last_update = time.time()
        self.frequency = frequency

    def update(self, dt):
        if time.time() - self.last_update > 1.0 / self.frequency:
            self.last_update = time.time()
            
            self.cpu_util = psutil.cpu_percent()
            self.cpu_util_history.append(self.cpu_util)
            self.cpu_util_history.pop(0)

            mem = psutil.virtual_memory()
            self.memory_util = mem.percent
            self.memory_available_mb = (mem.total - mem.available) / (1024 * 1024)
            self.memory_util_history.append(self.memory_util)
            self.memory_util_history.pop(0)
            self.fps = int(sum(self.fps_history)/len(self.fps_history))
        
        self.fps_history.append(1.0/dt)
        self.fps_history.pop(0)


    def draw(self, dt):
        self.update(dt)
        imgui.text(f"FPS: {self.fps}")
        imgui.text(f"CPU: {self.cpu_util}%")
        imgui.same_line()
        imgui.set_cursor_pos((imgui.get_window_width()*0.5+10, imgui.get_cursor_pos()[1]))
        imgui.text(f"Memory: {self.memory_available_mb:.1f}/{self.memory_total_mb:.1f}MB ({self.memory_util:.1f}%)")
        vals = array.array("f", self.cpu_util_history)
        imgui.plot_lines("", vals, scale_min=0.0, scale_max=100.0, graph_size=((imgui.get_window_width()-20)*0.5, 80))
        imgui.same_line()
        vals = array.array("f", self.memory_util_history)
        imgui.plot_lines("", vals, scale_min=0.0, scale_max=100.0, graph_size=((imgui.get_window_width()-20)*0.5, 80))
        

