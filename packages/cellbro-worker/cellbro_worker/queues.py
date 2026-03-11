def read_h5ad(file_path: str) -> str:
    from .tasks import io
    return io.read_h5ad.apply_async(args=[file_path])  # type: ignore

def qc(mt_prefix: str, ribo_prefixes: list[str], hb_pattern: str, percent_top: int) -> str:
    from .tasks import qc
    return qc.run.apply_async(args=[mt_prefix, ribo_prefixes, hb_pattern, percent_top])  # type: ignore

def test_plot(plot_id: str) -> str:
    from .tasks import test_plot
    return test_plot.apply_async(args=[plot_id])  # type: ignore


def plot_total_counts_histogram(plot_id: str) -> str:
    from .tasks import qc
    return qc.plot_total_counts_histogram.apply_async(args=[plot_id])  # type: ignore