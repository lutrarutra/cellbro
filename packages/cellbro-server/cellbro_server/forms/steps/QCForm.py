from fastapi import Response

from cellbro_worker import tasks

from ...core.cache import worker_redis
from ...core.HTMXForm import HTMXForm
from ...core.context import ctx
from ...components import inputs

class QCForm(HTMXForm):
    template_path = "forms/steps/qc.html"
    title = "Quality Control Parameters"

    mt_prefix = inputs.string.StringInputField("Mitochondrial Gene Prefix", placeholder="e.g., MT-", default="MT-")
    ribo_prefix = inputs.string.StringInputField("Ribosomal Gene Prefixes comma-separated", placeholder="e.g., RPS,RPL", default="RPS,RPL")
    hb_pattern = inputs.string.StringInputField("Hemoglobin Gene Pattern", placeholder="e.g., ^HB[^(P)]", default="^HB[^(P)]")

    percent_top = inputs.IntegerInputField("Percent Top Genes", placeholder="e.g., 50", default=20, min_value=0, max_value=100)

    async def process(self) -> Response:
        if not await self.validate():
            return await self.make_response()
        
        ribo_prefixes = [prefix.strip() for prefix in self.ribo_prefix.data.split(",") if prefix.strip()]
        
        task = await tasks.qc.run.kiq(
            mt_prefix=self.mt_prefix.data,
            ribo_prefixes=ribo_prefixes,
            hb_pattern=self.hb_pattern.data,
            percent_top=self.percent_top.data,
        )
        await worker_redis.client.set(f"task_redirect:{task.task_id}", ctx.url_for("qc_page"))
        return await self.loading_response(task_name="Quality Control", task_id=task.task_id)