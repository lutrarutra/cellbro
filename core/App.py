import multiprocessing

import imgui

import scanpy as sc

import util.FileType as FileType
import Window, util.Logger as Logger
import util.Dataset as Dataset
from analysis import preprocessing as pp
from interface import Dashboard
from util import Event

import analysis.projection.umap as projection_umap
import analysis.projection.trimap as projection_trimap
import analysis.projection.tsne as projection_tsne
import analysis.projection.pca as projection_pca
from analysis import violin

import io_window

class App():
    def __init__(self):
        self.running = True
        self.logger = Logger.Logger()
        self.window = Window.Window()
        self.dataset = None
        self.dashboard = Dashboard.Dashboard(self)
        self.event_handler = Event.EventHandler(self)
        self.processes = {}

    def loop(self):
        while True:
            if not self.window.prepare_frame():
                self.running = False
                return

            self.gui()
            self.window.render_frame()
            self.process_events()

    def process_events(self):
        for event_key in self.event_handler.completed_events:
            event = self.event_handler.events[event_key]
            if event_key == "load_file":
                self.dataset = Dataset.Dataset(path=event.args["path"], file_type=FileType.SC, logger=self.logger)
                self.dashboard.pipeline = pp.PreprocessProgress()
                self.event_handler.add_event("ask_annotate")
                self.dashboard.popup = pp.AskAnnotate(self.event_handler)
            
            elif event_key == "ask_annotate":
                if event.args["answer"]:
                    self.event_handler.add_event("load_annotation")
                    self.dashboard.popup = io_window.Io_window(self.event_handler, event_key="load_annotation")
                else:
                    self.event_handler.add_event("ask_preprocess")
                    self.dashboard.popup = pp.AskPreprocess(self.event_handler)
                    self.dashboard.pipeline.step += 1

            elif event_key == "load_annotation":
                self.dataset.load_annotation(path=event.args["path"])
                self.event_handler.add_event("ask_preprocess")
                self.dashboard.popup = pp.AskPreprocess(self.event_handler)
                self.dashboard.pipeline.step += 1

            elif event_key == "ask_preprocess":
                if event.args["answer"]:
                    self.event_handler.add_event("filter_form")
                    self.dashboard.main = pp.FilterForm(self.dataset, self.event_handler)

            elif event_key == "filter_form":
                self.event_handler.add_event("mtqc_form")
                self.dashboard.main = pp.MTQCForm(self.dataset, self.event_handler)
                self.dashboard.pipeline.step += 1
                self.processes["qc_mt_plot"] = multiprocessing.Process(target=sc.pl.scatter, args=(self.dataset.adata,), kwargs=dict(x="total_counts", y="pct_counts_mt", color="n_genes_by_counts"))
                self.processes["qc_mt_plot"].start()
            
            elif event_key == "mtqc_form":
                self.event_handler.add_event("normalize_form")
                self.dashboard.main = pp.NormalizeForm(self.dataset, self.event_handler)
                self.dashboard.pipeline.step += 1

            elif event_key == "normalize_form":
                self.dataset.post_process_init()
                self.dashboard.pipeline.step += 1

            self.event_handler.remove_event(event_key)

        self.event_handler.completed_events = []

    def gui(self):
        flags = 0
        if self.dashboard.popup is not None:
            if not self.dashboard.popup.draw():
                self.dashboard.popup = None
            else:
                flags |= imgui.WINDOW_NO_INPUTS
            
        if imgui.begin_main_menu_bar():
            if imgui.begin_menu("File", self.dashboard.popup==None):
                clicked_quit, _ = imgui.menu_item(
                    "Quit", 'Cmd+Q', False, True
                )

                clicked_load, _ = imgui.menu_item(
                    "Load", '', False, True
                )

                if clicked_quit:
                    exit(1)

                if clicked_load:
                    self.event_handler.add_event("load_file")
                    self.dashboard.popup = io_window.Io_window(self.event_handler, event_key="load_file")

                imgui.end_menu()

            if imgui.begin_menu("Projection", self.dataset is not None and self.dataset.preprocessed and self.dashboard.popup==None):

                if imgui.menu_item("UMAP", '', False, True)[0]:
                    self.dashboard.main = projection_umap.UMAP(self)

                if imgui.menu_item("Trimap", '', False, True)[0]:
                    self.dashboard.main = projection_trimap.Trimap(self)

                if imgui.menu_item("t-SNE", '', False, True)[0]:
                    self.dashboard.main = projection_tsne.TSNE(self)

                if imgui.menu_item("PCA", '', False, True)[0]:
                    self.dashboard.main = projection_pca.PCA(self)

                imgui.end_menu()

            if imgui.begin_menu("Plots", self.dataset is not None and self.dataset.preprocessed and self.dashboard.popup==None):
                if imgui.menu_item("Violin", '', False, True)[0]:
                    self.dashboard.main = violin.Violin(self)
                if imgui.menu_item("Heatmap", '', False, True)[0]:
                    pass
                imgui.end_menu()

            imgui.end_main_menu_bar()

        if self.dataset is not None:
            self.dashboard.pipeline.draw(flags=flags)

        # self.logger.draw(flags=flags)
        self.dashboard.draw_footer(flags=flags)

        imgui.begin("Main", flags=flags | imgui.WINDOW_NO_COLLAPSE | imgui.WINDOW_NO_TITLE_BAR)
        if self.dashboard.main:
            self.dashboard.main.draw()
        imgui.end()

        
