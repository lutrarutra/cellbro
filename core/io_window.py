import imgui
import os

import util.FileFormat as FileFormat
class Io_window():
    def __init__(self, event_handler, event_key="load_file", title="Select file", cwd=os.getcwd()):
        self.finished = False
        self.event_handler = event_handler
        self.event_key = event_key
        self.title = title
        self.path = ""
        self.show_hidden = False
        self.show_filetypes = 0
        self.cwd = ""
        self.set_cwd(cwd)

    def set_cwd(self, cwd):
        self.cwd = cwd
        self.paths = sorted([(f, os.path.isfile(os.path.join(self.cwd, f))) for f in os.listdir(cwd)], key=lambda x: x[0])

    def close(self):
        self.event_handler.remove_event(self.event_key)

    def draw(self):
        if self.finished:
            return False

        imgui.begin(self.title)

        imgui.text(self.cwd)
        imgui.begin_child("Files:", -10, -120, border=True)

        if imgui.button(".."):
            self.set_cwd(os.path.dirname(self.cwd))

        for path, is_file in self.paths:
            if path[0] == "." and not self.show_hidden:
                continue
            if is_file and self.show_filetypes == 0 and not "." + path.split(".")[-1] in FileFormat.sc_file_extensions:
                continue
            
            if is_file:
                if path[0] == ".":
                    imgui.push_style_color(imgui.COLOR_TEXT, 0.5, 0.5, 0.5)
                elif "." + path.split(".")[-1] in FileFormat.sc_file_extensions:
                    imgui.push_style_color(imgui.COLOR_TEXT, 0, 1, 0)
                else:
                    imgui.push_style_color(imgui.COLOR_TEXT, 1, 0.7, 0)
            else:
                imgui.push_style_color(imgui.COLOR_TEXT, 1, 1, 0)

            if imgui.button(f"{path}{'' if is_file else '/'} "):
                if not is_file:
                    self.set_cwd(os.path.join(self.cwd, path))
                else:
                    self.path = os.path.join(self.cwd, path)

            imgui.pop_style_color()
            

        imgui.end_child()

        _, self.show_filetypes = imgui.combo("Show file types", self.show_filetypes, ["/".join(FileFormat.sc_file_extensions), "*"])

        _, self.path = imgui.input_text(
            "File path", self.path, 256
        )
        _, self.show_hidden = imgui.checkbox("Show hidden", self.show_hidden)
        imgui.same_line()
        if self.path == "":
            imgui.push_style_color(imgui.COLOR_TEXT, 0.5, 0.5, 0.5)
        else:
            imgui.push_style_color(imgui.COLOR_TEXT, 1, 1, 1)
        if imgui.button("Open"):
            if self.path:
                imgui.pop_style_color()
                imgui.end()
                self.event_handler.complete_event(self.event_key, args=dict(path=self.path))
                return False
        imgui.pop_style_color()

        imgui.same_line()

        if imgui.button("Close"):
            imgui.end()
            self.event_handler.complete_event(self.event_key, args=dict(path=None))
            return False

        imgui.end()
        return True

