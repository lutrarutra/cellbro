import imgui
import os

import util.FileType as FileType

class Io_window():
    def __init__(self, title="Select file", cwd=os.getcwd()):
        self.opened = True
        self.title = title
        self.path = ""
        self.show_hidden = False
        self.show_filetypes = 0
        self.cwd = cwd
        self.paths = [(f, os.path.isfile(f)) for f in os.listdir(cwd)]

    def draw(self):
        if self.opened:
            imgui.begin(self.title)
            
            imgui.text(self.cwd)
            imgui.begin_child("Files:", -10, -70, border=True)
            if imgui.button(".."):
                self.cwd = os.path.dirname(self.cwd)
                self.paths = [(f, os.path.isfile(os.path.join(self.cwd, f))) for f in os.listdir(self.cwd)]

            for path, is_file in self.paths:
                if path[0] == "." and not self.show_hidden:
                    continue
                if is_file and self.show_filetypes == 0 and not path.split(".")[-1] in FileType.sc_file_extensions:
                    continue

                if imgui.button(f"{path}{'' if is_file else '/'} "):
                    if not is_file:
                        self.cwd = os.path.join(self.cwd, path)
                        self.paths = [(f, os.path.isfile(os.path.join(self.cwd, f))) for f in os.listdir(self.cwd)]
                    else:
                        self.path = os.path.join(self.cwd, path)

            imgui.end_child()

            clicked, self.show_filetypes = imgui.combo("Show file types", self.show_filetypes, ["/".join(FileType.sc_file_extensions), "*"])

            changed, text_val = imgui.input_text(
                "File path", self.path, 256
            )
            _, self.show_hidden = imgui.checkbox("Show hidden", self.show_hidden)
            imgui.same_line()
            imgui.dummy(imgui.get_window_width()-250, 0)
            imgui.same_line()
            imgui.button("Open", width=50)
            imgui.same_line()

            if imgui.button("Close", width=50):
                # self.opened = False
                imgui.end()
                return False

            imgui.end()
            return True

