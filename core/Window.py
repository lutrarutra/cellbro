import multiprocessing

import glfw
import OpenGL.GL as gl

import scanpy as sc

import imgui
from imgui.integrations.glfw import GlfwRenderer

import io_window
import util.FileType as FileType
from analysis import preprocessing as pp
import analysis.projection.umap as projection_umap
import analysis.projection.trimap as projection_trimap
import analysis.projection.tsne as projection_tsne
import analysis.projection.pca as projection_pca
from analysis import violin

class Window():
    def __init__(self, app, logger):
        imgui.create_context()
        imgui.get_io().font_global_scale = 2.0
        self.window = Window.impl_glfw_init()
        self.impl = GlfwRenderer(self.window)
        self.logger = logger
        self.childs = {}
        self.app = app

    def gui(self):
        if imgui.begin_main_menu_bar():
            if imgui.begin_menu("File", True):
                clicked_quit, selected_quit = imgui.menu_item(
                    "Quit", 'Cmd+Q', False, True
                )

                clicked_load, selected_load = imgui.menu_item(
                    "Load", '', False, True
                )

                if clicked_quit:
                    exit(1)

                if clicked_load:
                    if "file_browser" not in self.childs:
                        self.childs["file_browser"] = io_window.Io_window()

                imgui.end_menu()

            if imgui.begin_menu("Projection", self.app.dataset is not None and self.app.dataset.preprocessed):

                if imgui.menu_item("UMAP", '', False, True)[0]:
                    if "umap_form" not in self.childs.keys():
                        self.childs["umap_form"] = projection_umap.UMAP(self.app.dataset, self.app)

                if imgui.menu_item("Trimap", '', False, True)[0]:
                    if "trimap_form" not in self.childs.keys():
                        self.childs["trimap_form"] = projection_trimap.Trimap(self.app.dataset, self.app)

                if imgui.menu_item("t-SNE", '', False, True)[0]:
                    if "tsne_form" not in self.childs.keys():
                        self.childs["tsne_form"] = projection_tsne.TSNE(self.app.dataset, self.app)

                if imgui.menu_item("PCA", '', False, True)[0]:
                    if "pca_form" not in self.childs.keys():
                        self.childs["pca_form"] = projection_pca.PCA(self.app.dataset, self.app)

                imgui.end_menu()

            if imgui.begin_menu("Plots", self.app.dataset is not None and self.app.dataset.preprocessed):
                if imgui.menu_item("Violin", '', False, True)[0]:
                    if "violin_form" not in self.childs.keys():
                        self.childs["violin_form"] = violin.Violin(self.app.dataset, self.app)
                if imgui.menu_item("Heatmap", '', False, True)[0]:
                    pass
                imgui.end_menu()


            imgui.end_main_menu_bar()

        self.logger.draw()
        if self.app.dataset:
            self.app.dataset.draw()
        self.process_childs()
        
    def process_childs(self):
        for key in list(self.childs.keys()):
            opened, val = self.childs[key].draw()
            if not opened:
                del self.childs[key]

                if key == "file_browser" and val is not None:
                    self.logger.info(f"Loaded file: {val}")
                    self.app.init_dataset(val, FileType.SC)
                    self.logger.info(f"Cells: {self.app.dataset.adata.shape[0]} | Genes: {self.app.dataset.adata.shape[1]}")

                if key == "ask_phenodata":
                    if val:
                        self.childs["annotate_file_browser"] = io_window.Io_window()
                    else:
                        self.childs["ask_preprocessing"] = pp.AskPreprocessForm()

                
                if key == "ask_preprocessing":
                    if val:
                        self.childs["filter_form"] = pp.FilterForm(self.app.dataset)
                        self.childs["preprocessing_progress"].step += 1
                    else:
                        self.childs["preprocessing_progress"].finished = True
                        self.app.dataset.post_process_init()

                if key == "annotate_file_browser" and val is not None:
                    self.logger.info(f"Loaded file: {val}")
                    self.app.dataset.annotate(val)
                    self.childs["filter_form"] = pp.FilterForm(self.app.dataset)
                    self.childs["preprocessing_progress"].step += 1


                if key == "filter_form":
                    self.logger.info(f"Dataset filtered...")
                    self.logger.info(f"Cells: {self.app.dataset.adata.shape[0]} | Genes: {self.app.dataset.adata.shape[1]}")
                    if "qc_mt_plot" not in self.app.processes.keys():
                        self.app.processes["qc_mt_plot"] = multiprocessing.Process(target=sc.pl.scatter, args=(self.app.dataset.adata,), kwargs=dict(x="total_counts", y="pct_counts_mt", color="n_genes_by_counts"))
                        self.app.processes["qc_mt_plot"].start()
                        self.childs["qc_mt_form"] = pp.MTQCForm(self.app.dataset)
                        self.childs["preprocessing_progress"].step += 1

                if key == "qc_mt_form":
                    self.logger.info("MT QC finished...")
                    self.logger.info(f"Cells: {self.app.dataset.adata.shape[0]} | Genes: {self.app.dataset.adata.shape[1]}")
                    self.childs["normalize_form"] = pp.NormalizeForm(self.app.dataset)
                    self.childs["preprocessing_progress"].step += 1


                if key == "normalize_form":
                    self.logger.info("Created layers: X (log1p), 'counts', 'ncounts', 'centered', 'logcentered', neighbors and 'X_pca'")
                    self.logger.info("Normalization finished...")
                    self.childs["preprocessing_progress"].finished = True
                    self.app.dataset.post_process_init()

                if key == "umap_form":
                    if val is None:
                        continue
                    val.apply()
                
                if key == "trimap_form":
                    if val is None:
                        continue
                    val.apply()

                if key == "tsne_form":
                    if val is None:
                        continue
                    val.apply()
                
                if key == "pca_form":
                    if val is None:
                        continue
                    val.apply()
                

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

        gl.glClearColor(0.0,0.0,0.0,1.0)
        gl.glClear(gl.GL_COLOR_BUFFER_BIT)

        imgui.render()
        self.impl.render(imgui.get_draw_data())
        glfw.swap_buffers(self.window)
        return True

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