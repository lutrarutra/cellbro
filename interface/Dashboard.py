import imgui

import math

class Spinner():
    def __init__(self, window):
        self.window = window
        self.stage = 0
        self.velocity = 7

    def draw(self, dt):
        self.stage += dt * self.velocity
        w, h = self.window.get_window_size()
        x = w * 0.5
        y = h * 0.5-20
        a3 = math.cos(self.stage) * 0.5 + 0.5
        a2 = math.cos(self.stage + math.pi / 3) * 0.5 + 0.5
        a1 = math.cos(self.stage + 2 * math.pi / 3) * 0.5 + 0.5
        draw_list = imgui.get_window_draw_list()
        draw_list.add_circle_filled(x-25,   y, 10, imgui.get_color_u32_rgba(1.0, 1.0, 1.0, a1), 32)
        draw_list.add_circle_filled(x,      y, 10, imgui.get_color_u32_rgba(1.0, 1.0, 1.0, a2), 32)
        draw_list.add_circle_filled(x+25,   y, 10, imgui.get_color_u32_rgba(1.0, 1.0, 1.0, a3), 32)

class LoadingPopup():
    def __init__(self, window):
        self.window = window
        self.spinner = Spinner(self.window)
        self.w = 400
        self.h = 200

    def draw(self, dt, flags=imgui.WINDOW_NO_COLLAPSE | imgui.WINDOW_NO_RESIZE | imgui.WINDOW_NO_TITLE_BAR):
        w, h = self.window.get_window_size()
        imgui.set_next_window_position((w - self.w) * 0.5 , (h - self.h) * 0.5)
        imgui.set_next_window_size(self.w, self.h)
        imgui.begin("", flags=flags)
        self.spinner.draw(dt)
        w,h = imgui.get_window_size()
        text_w, text_h = imgui.calc_text_size("Loading...")
        imgui.set_cursor_pos((w*0.5 - text_w*0.5, h*0.5 - text_h*0.5+20))
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

    def draw_footer(self, dt, flags=0):
        imgui.begin("Footer", flags=flags | imgui.WINDOW_NO_TITLE_BAR)

        imgui.columns(2 + (1 if self.app.dataset else 0), "tabs")
        if imgui.selectable("Log", self.footer_tab == 0)[0]:
            self.footer_tab = 0
        imgui.next_column()
        if imgui.selectable(f"Performance", self.footer_tab == 1)[0]:
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
            self.app.htop.draw(dt=dt)
        elif self.footer_tab == 2:
            if self.app.dataset:
                self.app.dataset.draw()

        imgui.end()