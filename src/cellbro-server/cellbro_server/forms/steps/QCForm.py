from fastapi import Response

from cellbro_worker import queues

from ...core.responses import htmx_response
from ...core.HTMXForm import HTMXForm
from ...core.context import ctx
from ...components import inputs

class QCForm(HTMXForm):
    template_path = "forms/steps/qc.html"

    mt_prefix = inputs.string.StringInputField("Mitochondrial Gene Prefix", placeholder="e.g., MT-", default="MT-")
    ribo_prefix = inputs.string.StringInputField("Ribosomal Gene Prefixes comma-separated", placeholder="e.g., RPS,RPL", default="RPS,RPL")
    hb_pattern = inputs.string.StringInputField("Hemoglobin Gene Pattern", placeholder="e.g., ^HB[^(P)]", default="^HB[^(P)]")

    percent_top = inputs.IntegerInputField("Percent Top Genes", placeholder="e.g., 50", default=20, min_value=0, max_value=100)

    async def process(self) -> Response:
        if not await self.validate():
            return await self.make_response()
        
        ribo_prefixes = [prefix.strip() for prefix in self.ribo_prefix.data.split(",") if prefix.strip()]
        
        queues.qc(
            mt_prefix=self.mt_prefix.data,
            ribo_prefixes=ribo_prefixes,
            hb_pattern=self.hb_pattern.data,
            percent_top=self.percent_top.data,
        )
        return await htmx_response(redirect=ctx.request.url_for("qc_page"))