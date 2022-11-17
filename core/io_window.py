import imgui
import os

import util.FileFormat as FileFormat
class Io_window():
    def __init__(self, title="Select file", cwd=os.getcwd()):
        self.opened = True
        self.title = title
        self.path = ""
        self.show_hidden = False
        self.show_filetypes = 0
        self.cwd = ""
        self.set_cwd(cwd)

    def set_cwd(self, cwd):
        self.cwd = cwd
        self.paths = sorted([(f, os.path.isfile(os.path.join(self.cwd, f))) for f in os.listdir(cwd)], key=lambda x: x[0])

    def draw(self):
        if self.opened:
            imgui.begin(self.title)
            
            imgui.text(self.cwd)
            imgui.begin_child("Files:", -10, -70, border=True)

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

            clicked, self.show_filetypes = imgui.combo("Show file types", self.show_filetypes, ["/".join(FileFormat.sc_file_extensions), "*"])

            changed, text_val = imgui.input_text(
                "File path", self.path, 256
            )
            _, self.show_hidden = imgui.checkbox("Show hidden", self.show_hidden)
            imgui.same_line()
            imgui.dummy(imgui.get_window_width()-250, 0)
            imgui.same_line()

            if imgui.button("Open", width=50):
                imgui.end()
                return False, self.path

            imgui.same_line()

            if imgui.button("Close", width=50):
                # self.opened = False
                imgui.end()
                return False, None

            imgui.end()
            return True, None

