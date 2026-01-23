from fastapi import Depends, APIRouter
import sqlalchemy as sa

router = APIRouter(prefix="/htmx/steps", tags=["steps", "htmx"])

from cellbro_db import models
from cellbro_db.core.session import AsyncSession

from ...core.context import ctx
from ...core import responses
from ... import forms, logic
from ...core.dependencies import db_session, dataset


@router.get("/qc")
@router.post("/qc")
async def qc_step(db: AsyncSession = Depends(db_session), _ = Depends(dataset)):
    form = forms.steps.QCForm()

    if ctx.request.method == "POST":
        return await form.process()
    
    return await form.make_response()
    