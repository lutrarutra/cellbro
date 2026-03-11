from fastapi import Response

from cellbro_db.core.session import AsyncSession

from ... import logic
from ...core import exceptions as exc
from ...core.responses import htmx_response
from ...core.HTMXForm import HTMXForm
from ...core.context import ctx
from ...components import inputs

class LoginForm(HTMXForm):
    template_path = "forms/auth/login.html"

    email = inputs.string.StringInputField("Email", placeholder="Enter your email")
    password = inputs.string.PasswordInputField("Password",  placeholder="Enter your password")
    
    async def validate(self, db: AsyncSession) -> bool:
        if not await super().validate():
            return False
        
        try:
            await logic.auth.login(
                db=db, 
                email=self.email.data, 
                password=self.password.data
            )
        except exc.InvalidCredentialsException:
            self.email.errors.append("Invalid email or password")
            self.password.errors.append("Invalid email or password")
            return False

        return True
    
    async def process(self, db: AsyncSession) -> Response:
        if not await self.validate(db):
            return await self.make_response()
        
        user = await logic.auth.login(db=db, email=self.email.data, password=self.password.data)
        await ctx.flash("Login Successful!", category="success")
        return await htmx_response(redirect="/")