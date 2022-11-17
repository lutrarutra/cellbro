import glfw
import OpenGL.GL as gl

import imgui
from imgui.integrations.glfw import GlfwRenderer

import io_window

class Window():
    def __init__(self):
        imgui.create_context()
        self.window = Window.impl_glfw_init()
        self.impl = GlfwRenderer(self.window)

        self.childs = {}

    def gui(self):
        if imgui.begin_main_menu_bar():
            if imgui.begin_menu("File", True):

                clicked_quit, selected_quit = imgui.menu_item(
                    "Quit", 'Cmd+Q', False, True
                )

                clicked_load, selected_load = imgui.menu_item(
                    "Load", 'Cmd+L', False, True
                )

                if clicked_quit:
                    exit(1)

                if clicked_load:
                    if "file_browser" not in self.childs:
                        self.childs["file_browser"] = io_window.Io_window()

                imgui.end_menu()

            imgui.end_main_menu_bar()
        
        for key in list(self.childs.keys()):
            if not self.childs[key].draw():
                del self.childs[key]
    
    def update(self):
        if glfw.window_should_close(self.window):
            self.impl.shutdown()
            glfw.terminate()
            return False

        glfw.poll_events()
        self.impl.process_inputs()

        imgui.new_frame()
        self.gui()
        
        #imgui.show_test_window()

        gl.glClearColor(0.3,0.3,0.3,0.1)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT)

        imgui.render()
        self.impl.render(imgui.get_draw_data())
        glfw.swap_buffers(self.window)
        return True

    @staticmethod 
    def impl_glfw_init():
        width, height = 1280, 720
        window_name = "CellBro"

        if not glfw.init():
            print("Could not initialize OpenGL context")
            exit(1)

        # OS X supports only forward-compatible core profiles from 3.2
        glfw.window_hint(glfw.CONTEXT_VERSION_MAJOR, 3)
        glfw.window_hint(glfw.CONTEXT_VERSION_MINOR, 3)
        glfw.window_hint(glfw.OPENGL_PROFILE, glfw.OPENGL_CORE_PROFILE)
        glfw.window_hint(glfw.RESIZABLE, 1)
        glfw.window_hint(glfw.TRANSPARENT_FRAMEBUFFER, 1)

        glfw.window_hint(glfw.OPENGL_FORWARD_COMPAT, gl.GL_TRUE)

        # Create a windowed mode window and its OpenGL context
        window = glfw.create_window(
            int(width), int(height), window_name, None, None
        )
        glfw.make_context_current(window)

        if not window:
            glfw.terminate()
            print("Could not initialize Window")
            exit(1)

        return window