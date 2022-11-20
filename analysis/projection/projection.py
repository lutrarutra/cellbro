import webbrowser

import imgui

from util import Query
class Projection():
    def __init__(self, app):
        self.app = app
        self.type = "Projection"
        self.finished = False
        self.leiden = True
        self.current_tab = 0

        self.key_query = Query.Query(
            sorted(list(self.app.dataset.adata.obs.columns)) + sorted(list(self.app.dataset.adata.var.index)),
            proposal_keys=list(self.app.dataset.adata.obs.columns)
        )

        self.selected_keys = []
        self.calc_params = {}
        self.plot_params = {
            "ncols": 1,
            "frameon": False,
        }

        self.leiden_params = {
            "resolution": 1.0,
            "random_state": 0,
        }

        self.query = ""

    def draw(self):
        if self.finished:
            return False

        imgui.columns(4, "projection")
        if imgui.selectable(self.type, self.current_tab == 0)[0]:
            self.current_tab = 0
        imgui.next_column()
        if imgui.selectable("Leiden", self.current_tab == 1)[0]:
            self.current_tab = 1
        imgui.next_column()
        if imgui.selectable("Features", self.current_tab == 2)[0]:
            self.current_tab = 2
        imgui.next_column()
        if imgui.selectable("Plotting", self.current_tab == 3)[0]:
            self.current_tab = 3
        imgui.columns(1)
        imgui.separator()

        if self.current_tab == 0:
            for key, value in self.calc_params.items():
                if isinstance(value, int):
                    _, self.calc_params[key] = imgui.input_int(key, value)
                elif isinstance(value, float):
                    _, self.calc_params[key] = imgui.input_float(key, value)

            imgui.set_cursor_pos((imgui.get_cursor_pos()[0], imgui.get_window_height() - 40))
            if imgui.button("Next"):
                self.current_tab += 1

        elif self.current_tab == 1:
            imgui.same_line()
            _, self.leiden = imgui.checkbox("Enabled", self.leiden)
            if self.leiden:
                for key, value in self.leiden_params.items():
                    if isinstance(value, int):
                        _, self.leiden_params[key] = imgui.input_int(key, value)
                    elif isinstance(value, float):
                        _, self.leiden_params[key] = imgui.input_float(key, value)
            
            imgui.set_cursor_pos((imgui.get_cursor_pos()[0], imgui.get_window_height() - 40))
            if imgui.button("Next"):
                self.current_tab += 1
        
        elif self.current_tab == 2:
            imgui.text("Search:")
            imgui.push_item_width(imgui.get_window_width()*0.5)
            query_changed, self.query  = imgui.input_text(
                    "Key (Feature/GeneName)", self.query, 256
            )
            imgui.pop_item_width()
            if query_changed:
                self.key_query.query(self.query)

            imgui.begin_child("Available features text", (imgui.get_window_width()-30)*0.5, 40, border=False)
            imgui.text("Available Features:")
            imgui.end_child()
            imgui.same_line()
            imgui.begin_child("Selected features text", (imgui.get_window_width()-30)*0.5, 40, border=False)
            imgui.text("Selected Features:")
            imgui.end_child()
            imgui.begin_child("Available Features", (imgui.get_window_width()-30)*0.5, -50, border=True)
            for key in self.key_query.proposal_keys:
                clicked, _ = imgui.selectable(key, key in self.selected_keys)
                if clicked:
                    if key in self.selected_keys:
                        self.selected_keys.remove(key)
                    else:
                        self.selected_keys.append(key)

            imgui.end_child()
            imgui.same_line()
            imgui.begin_child("Selected Features", (imgui.get_window_width()-30)*0.5, -50, border=True)
            for key in self.selected_keys:
                if imgui.button(key):
                    self.selected_keys.remove(key)

            imgui.end_child()
            imgui.set_cursor_pos((imgui.get_cursor_pos()[0], imgui.get_window_height() - 40))
            if imgui.button("Next"):
                self.current_tab += 1

        elif self.current_tab == 3:
            for key, value in self.plot_params.items():
                if isinstance(value, int):
                    _, self.plot_params[key] = imgui.input_int(key, value)
                elif isinstance(value, float):
                    _, self.plot_params[key] = imgui.input_float(key, value)
                elif isinstance(value, bool):
                    _, self.plot_params[key] = imgui.checkbox(key, value)

            if self.type == "PCA":
                self.ask_dimensions()

            imgui.set_cursor_pos((imgui.get_cursor_pos()[0], imgui.get_window_height() - 40))
            if imgui.button("Plot"):
                self.finished = True
                self.apply()
                return False

        imgui.same_line()
        if imgui.button("Cancel"):
            return False

        if self.current_tab == 0:
            imgui.same_line()
            if imgui.button("Documentation"):
                webbrowser.open(self.tl_documentation)
        elif self.current_tab == 1:
            imgui.same_line()
            if imgui.button("Documentation"):
                webbrowser.open("https://scanpy.readthedocs.io/en/stable/generated/scanpy.tl.leiden.html")
        elif self.current_tab == 3:
            imgui.same_line()
            if imgui.button("Documentation"):
                webbrowser.open(self.pl_documentation)

        return True