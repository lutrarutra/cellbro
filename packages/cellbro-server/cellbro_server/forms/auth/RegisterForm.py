from fastapi import Response

from cellbro_db.core.session import AsyncSession
from cellbro_db import types, queries

from ... import logic
from ...core import exceptions as exc
from ...core.responses import htmx_response
from ...core.HTMXForm import HTMXForm
from ...core.context import ctx
from ...components import inputs

class RegisterForm(HTMXForm):
    template_path = "forms/auth/register.html"

    first_name = inputs.string.StringInputField(
        "First Name", placeholder="Enter your first name",
        min_length=1
    )
    last_name = inputs.string.StringInputField(
        "Last Name", placeholder="Enter your last name",
        min_length=1
    )

    email = inputs.string.EmailInputField(
        "Email", placeholder="Enter your email",
        autocomplete="username",
    )
    password = inputs.string.PasswordInputField(
        "Password",  placeholder="Enter your password",
        min_length=8, autocomplete="new-password"
    )
    repeat_password = inputs.string.PasswordInputField(
        "Repeat Password", placeholder="Repeat your password",
        min_length=8, autocomplete="new-password"
    )
    
    async def validate(self, db: AsyncSession) -> bool:
        if not await super().validate():
            return False
        
        if self.password.data != self.repeat_password.data:
            self.password.errors.append("Passwords do not match")
            self.repeat_password.errors.append("Passwords do not match")
            return False

        if (user := await db.get_one(queries.user.get(email=self.email.data))) is not None:
            self.email.errors.append("Email is already registered")
            return False
        
        return True
    
    async def process(self, db: AsyncSession) -> Response:
        if not await self.validate(db):
            return await self.make_response()
        
        user = await logic.auth.register(
            db=db, email=self.email.data, password=self.password.data,
            first_name=self.first_name.data, last_name=self.last_name.data,
            account_type=types.UserType.REGULAR
        )
        await ctx.flash("Registration Successful!", category="success")
        return await htmx_response(redirect="/")