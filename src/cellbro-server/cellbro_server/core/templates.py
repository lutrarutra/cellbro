import jinjax
from datetime import datetime
from fastapi.templating import Jinja2Templates
import markdown
from markupsafe import Markup
import pandas as pd

from cellbro_db import types

from .context import ctx
from .config import settings

templates = Jinja2Templates(directory="templates")

templates.env.add_extension(jinjax.JinjaX)
catalog = jinjax.Catalog(jinja_env=templates.env)
catalog.add_folder("templates/components")

templates.env.globals["types"] = types

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

def format_number(value: int | float | None) -> str:
    if pd.isna(value):
        return ""
    
    return f"{value:,}".replace(",", " ")

templates.env.filters["format_timestamp"] = format_timestamp
templates.env.filters["format_datetime"] = format_datetime
templates.env.filters["format_markdown"] = format_markdown
templates.env.filters["filesize"] = filesize
templates.env.filters["format_number"] = format_number


def render_component(component: str, **context) -> str:
    return catalog.render(component, **context)