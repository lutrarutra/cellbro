from celery.result import AsyncResult

from cellbro_db import types

def read_h5ad(file_path: str) -> AsyncResult:
    from .tasks import read_h5ad
    return read_h5ad.apply_async(args=[file_path])  # type: ignore

def qc(mt_prefix: str, ribo_prefixes: list[str], hb_pattern: str, percent_top: int) -> AsyncResult:
    from .tasks import qc
    return qc.apply_async(args=[mt_prefix, ribo_prefixes, hb_pattern, percent_top])  # type: ignore


def test_plot(plot_id: str) -> AsyncResult:
    from .tasks import test_plot
    return test_plot.apply_async(args=[plot_id])  # type: ignore


def plot_total_counts_histogram(plot_id: str) -> AsyncResult:
    from .tasks import plot_total_counts_histogram
    return plot_total_counts_histogram.apply_async(args=[plot_id])  # type: ignore