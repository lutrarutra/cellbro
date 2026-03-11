import json

from ..broker import task_broker
from . import io, qc

@task_broker.task
def test_plot(plot_id: str):
    from bokeh.plotting import figure
    from bokeh.themes import built_in_themes
    from bokeh.document import Document
    from bokeh.embed import json_item

    import time
    time.sleep(1)

    p = figure(x_axis_label="X", y_axis_label="Y", sizing_mode="stretch_both", match_aspect=True, output_backend="webgl")
    p.line([1, 2, 3, 4, 5], [6, 7, 2, 4, 5], line_width=2, color="navy")

    doc = Document()
    doc.add_root(p)
    doc.theme = built_in_themes['dark_minimal']

    self.r.publish(f'plot:{plot_id}', json.dumps(json_item(p, plot_id)))  # type: ignore


