from fastapi import Depends, APIRouter

from cellbro_db import models
from cellbro_db.core.session import AsyncSession

from ... import forms, logic
from ...core.context import ctx
from ...core.responses import htmx_response
from ...core.dependencies import db_session, get_user

router = APIRouter(prefix="/htmx/auth", tags=["auth", "htmx"])

@router.get("/login")
@router.post("/login")
async def htmx_login(db: AsyncSession = Depends(db_session)):
    form = forms.auth.LoginForm()
    if ctx.request.method == "POST":
        return await form.process(db)    
    return await form.make_response()


@router.post("/logout")
async def logout(user: models.User = Depends(get_user)):
    await logic.auth.logout(user.id)
    await ctx.flash("Logged Out!", category="info")
    return await htmx_response(redirect="/auth/login")


@router.get("/register")
@router.post("/register")
async def register(db: AsyncSession = Depends(db_session)):
    form = forms.auth.RegisterForm()
    if ctx.request.method == "POST":
        return await form.process(db)
    return await form.make_response()