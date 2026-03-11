from fastapi import Response

from cellbro_worker import queues
from cellbro_db import types

from ...core.responses import htmx_response
from ...core.HTMXForm import HTMXForm
from ...core.context import ctx
from ...components import inputs

class ScatterFigureForm(HTMXForm):
    template_path = "forms/figures/scatter.html"
    
    level = inputs.SelectInput("Level", options=[(t.value, t.name.title()) for t in types.AnnDataLayerType])

    
