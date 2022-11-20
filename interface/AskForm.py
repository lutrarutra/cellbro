import imgui

from interface import Form

import util.gui

class AskForm(Form.Form):
    def __init__(self, event_handler, event_key, question, question_tip=None, options={"Yes":True, "No":False}):
        super().__init__(event_handler, event_key)
        self.question = question
        self.options = options
        self.question_tip = question_tip

    def draw(self, flags=imgui.WINDOW_NO_MOVE | imgui.WINDOW_NO_COLLAPSE):
        if self.finished:
            return False
        imgui.begin("Question", flags=flags)
        w,h = imgui.get_window_size()
        imgui.set_cursor_position((w*0.5 - imgui.calc_text_size(self.question)[0]*0.5, h*0.5 - imgui.calc_text_size(self.question)[1]*0.5))
        imgui.text(self.question)
        if self.question_tip is not None:
            imgui.same_line()
            util.gui.tooltip(self.question_tip)

        for i, (key, val) in enumerate(self.options.items()):
            if type(key) == str:
                imgui.set_cursor_position((w*0.5 - len(self.options)*70 + (i) * 140 + 35, h-50))
                if imgui.button(key, 70, 30):
                    self.event_handler.complete_event(self.event_key, dict(answer=val))
                    self.finished = True
                    imgui.end()
                    return False
            else:
                # TODO: Add support for other types of options
                assert False

        imgui.end()
        return True