import jinjax
from datetime import datetime
from fastapi.templating import Jinja2Templates
import markdown
from markupsafe import Markup

from .context import ctx
from .config import settings

templates = Jinja2Templates(directory="templates")

templates.env.add_extension(jinjax.JinjaX)
catalog = jinjax.Catalog(jinja_env=templates.env)
catalog.add_folder("templates/components")

# Configure filters
def format_timestamp(value: float, format: str = "%Y-%m-%d %H:%M"):
    """Convert timestamp to readable string"""
    return datetime.fromtimestamp(value).strftime(format)

def format_datetime(value: datetime, format: str = "%Y-%m-%d %H:%M"):
    """Convert datetime to readable string"""
    return value.astimezone(settings.TIMEZONE).strftime(format)

def filesize(value: int | float):
    """Convert bytes to human readable format"""
    for unit in ['B', 'KB', 'MB', 'GB', 'TB']:
        if value < 1024.0:
            if unit == 'B':
                return f"{int(value)} {unit}"
            return f"{value:.1f} {unit}"
        value /= 1024.0
    return f"{value:.1f} PB"

def format_markdown(value: str) -> str:
    return Markup(markdown.markdown(value))

catalog.jinja_env.filters["format_timestamp"] = format_timestamp
catalog.jinja_env.filters["format_datetime"] = format_datetime
catalog.jinja_env.filters["format_markdown"] = format_markdown
catalog.jinja_env.filters["filesize"] = filesize


def render_component(component: str, **context) -> str:
    return catalog.render(component, **context)