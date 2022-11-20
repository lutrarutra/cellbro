import imgui, glfw

class LoadingPopup():
    def __init__(self, window):
        self.window = window
        self.w = 400
        self.h = 200

    def draw(self, flags=imgui.WINDOW_NO_COLLAPSE | imgui.WINDOW_NO_RESIZE | imgui.WINDOW_NO_TITLE_BAR):
        w, h = self.window.get_window_size()
        imgui.set_next_window_position((w - self.w) * 0.5 , (h - self.h) * 0.5)
        imgui.set_next_window_size(self.w, self.h)
        imgui.begin("", flags=flags)
        imgui.set_cursor_pos((imgui.get_window_width() * 0.5 - 60, imgui.get_window_height() * 0.5-10))
        imgui.text("Loading...")
        imgui.end()

class Dashboard():
    def __init__(self, app):
        self.app = app
        self.main = None
        self.pipeline = None
        self.popup = None
        self.footer = None
        self.footer_tab = 0
        self.blocking_popup = LoadingPopup(self.app.window)

    def draw_footer(self, flags=0):
        imgui.begin("Footer", flags=flags | imgui.WINDOW_NO_TITLE_BAR)

        imgui.columns(2 + (1 if self.app.dataset else 0), "tabs")
        if imgui.selectable("Log", self.footer_tab == 0)[0]:
            self.footer_tab = 0
        imgui.next_column()
        if imgui.selectable("Performance", self.footer_tab == 1)[0]:
            self.footer_tab = 1
        imgui.next_column()
        if self.app.dataset:
            if imgui.selectable("Dataset", self.footer_tab == 2)[0]:
                self.footer_tab = 2
            # imgui.next_column()
        imgui.columns(1)
        imgui.separator()

        if self.footer_tab == 0:
            self.app.logger.draw()
        elif self.footer_tab == 1:
            self.app.htop.draw()
        elif self.footer_tab == 2:
            if self.app.dataset:
                self.app.dataset.draw()

        imgui.end()