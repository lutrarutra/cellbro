import imgui

def tooltip(text, same_line=True):
    if same_line:
        imgui.same_line()
    imgui.text("(?)")
    if imgui.is_item_hovered():
        imgui.begin_tooltip()
        imgui.text(text)
        imgui.end_tooltip()