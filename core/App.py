import multiprocessing, threading, time

import imgui

import scanpy as sc

import util.htop as htop
import util.FileType as FileType
import Window, util.Logger as Logger
import util.Dataset as Dataset
from analysis import preprocessing as pp
from interface import Dashboard
from util import Event
from util import Task

import analysis.projection.UMAP as projection_umap
import analysis.projection.Trimap as projection_trimap
import analysis.projection.TSNE as projection_tsne
import analysis.projection.PCA as projection_pca
from analysis import Violin
from analysis import Heatmap

import io_window

class App():
    def __init__(self):
        self.running = True
        self.logger = Logger.Logger()
        self.window = Window.Window()
        self.dataset = None
        self.dashboard = Dashboard.Dashboard(self)
        self.event_handler = Event.EventHandler(self)
        self.task_handler = Task.TaskHandler()
        self.htop = htop.HTOP()
        self.figures = {}
        self.last_frame = time.time()
        self.dt = 0.0

    def loop(self):
        while True:
            self.dt = time.time() - self.last_frame
            self.last_frame = time.time()
            if not self.window.prepare_frame():
                self.running = False
                return

            self.gui()
            self.window.render_frame()
            self.process_events()
            self.task_handler.process_tasks()

    def process_events(self):
        for event_key in self.event_handler.completed_events:
            event = self.event_handler.events[event_key]
            if event_key == "load_file":
                self.dataset = Dataset.Dataset(path=event.args["path"], file_type=FileType.SC, logger=self.logger)
                self.task_handler.add_task("load_file", Task.Task(target=self.dataset.load_file))
                self.task_handler.tasks["load_file"].start()

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
                if event.args["path"] is not None:
                    self.task_handler.add_task("load_annotation", Task.Task(target=self.dataset.load_annotation, args=(event.args["path"],)))
                    self.task_handler.tasks["load_annotation"].start()

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

                self.figures["qc_mt_plot"] = multiprocessing.Process(
                    target=sc.pl.scatter, kwargs=dict(adata=self.dataset.adata, x="total_counts", y="pct_counts_mt", color="n_genes_by_counts")
                )
                self.figures["qc_mt_plot"].start()
            
            elif event_key == "mtqc_form":
                self.event_handler.add_event("normalize_form")
                self.dashboard.main = pp.NormalizeForm(self.dataset, self.event_handler)
                self.dashboard.pipeline.step += 1

            elif event_key == "normalize_form":
                self.task_handler.add_task("post_process_init", Task.Task(target=self.dataset.post_process_init))
                self.task_handler.tasks["post_process_init"].start()
                self.dashboard.pipeline.step += 1

            self.event_handler.remove_event(event_key)

        self.event_handler.completed_events = []

    def gui(self):
        flags = 0
        if self.task_handler.is_blocking():
            self.dashboard.blocking_popup.draw(dt=self.dt)
            flags |= imgui.WINDOW_NO_INPUTS
            

        elif self.dashboard.popup is not None:
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
                    self.dashboard.main = Violin.Violin(self)
                if imgui.menu_item("Heatmap", '', False, True)[0]:
                    self.dashboard.main = Heatmap.Heatmap(self)
                imgui.end_menu()

            imgui.end_main_menu_bar()

        if self.dataset is not None:
            self.dashboard.pipeline.draw(flags=flags)

        self.dashboard.draw_footer(flags=0, dt=self.dt)

        imgui.begin("Main", flags=flags | imgui.WINDOW_NO_COLLAPSE | imgui.WINDOW_NO_TITLE_BAR)
        if self.dashboard.main:
            if not self.dashboard.main.draw():
                self.dashboard.main = None
        imgui.end()

        
