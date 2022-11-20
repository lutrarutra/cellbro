import imgui

class Dashboard():
    def __init__(self, app):
        self.app = app
        self.main = None
        self.pipeline = None
        self.popup = None
        self.footer = None
        self.footer_tab = 0

    def draw_footer(self, flags=0):
        imgui.begin("Footer", flags=flags | imgui.WINDOW_NO_TITLE_BAR)

        imgui.columns(1 + (1 if self.app.dataset else 0), "tabs")
        if imgui.selectable("Log", self.footer_tab == 0)[0]:
            self.footer_tab = 0
        imgui.next_column()
        if self.app.dataset:
            if imgui.selectable("Dataset", self.footer_tab == 1)[0]:
                self.footer_tab = 1
            # imgui.next_column()
        imgui.columns(1)
        imgui.separator()

        if self.footer_tab == 0:
            self.app.logger.draw()
        elif self.footer_tab == 1:
            if self.app.dataset:
                self.app.dataset.draw()

        imgui.end()