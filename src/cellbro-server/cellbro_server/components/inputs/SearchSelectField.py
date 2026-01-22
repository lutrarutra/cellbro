from starlette.datastructures import URL

from .InputField import InputField

from .InputContext import InputContext

class SearchSelectField(InputContext):
    def __init__(
        self,
        label: str,
        search_url: URL,
        placeholder: str | None = "Search",
        default: str | None = None,
        name: str | None = None,
        required: bool = True,
    ):
        super().__init__(component="search.SearchSelect")
        self.selected_id = InputField(
            name=name or label.lower().replace(" ", "_") + "_id",
            label=label,
            default=default,
            pydantic_type=int,
            type="hidden",
            required=required,
        )
        self.selected_name = InputField(
            name=name or label.lower().replace(" ", "_") + "_name",
            label=label,
            default=default,
            pydantic_type=str,
            type="hidden",
            required=required,
        )
        self.placeholder = placeholder
        self.search_url = search_url