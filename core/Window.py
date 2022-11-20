import glfw
import OpenGL.GL as gl

import imgui
from imgui.integrations.glfw import GlfwRenderer

class Window():
    def __init__(self):
        imgui.create_context()
        imgui.get_io().font_global_scale = 2.0
        self.window = Window.impl_glfw_init()
        self.impl = GlfwRenderer(self.window)

    def prepare_frame(self):
        if glfw.window_should_close(self.window):
            self.impl.shutdown()
            glfw.terminate()
            return False

        glfw.poll_events()
        self.impl.process_inputs()
        imgui.new_frame()
        return True
    
    def render_frame(self):
        gl.glClearColor(0.0,0.0,0.0,1.0)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT)

        imgui.render()
        self.impl.render(imgui.get_draw_data())
        glfw.swap_buffers(self.window)

    @staticmethod 
    def impl_glfw_init():
        # width, height = 1280, 720
        width, height = 1920, 1080
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

    def get_window_size(self):
        return glfw.get_framebuffer_size(self.window)