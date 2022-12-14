from dataclasses import dataclass
class ParamsDict():
    def __init__(self, params_list):
        self.params = {}
        for param in params_list:
            self.params[param.key] = param

    def update(self, params_dict):
        for key in params_dict.keys():
            self.params[key].value = params_dict[key]

        return self

    def __getitem__(self, key):
        return self.params[key]

    def add(self, param):
        self.params[param.key] = param

    def pop(self, key):
        return self.params.pop(key)

    def remove(self, key):
        del self.params[key]

    def keys(self):
        return self.params.keys()

    def values(self):
        return self.params.values()

    def items(self):
        return self.params.items()

    def unravel(self):
        params = {}
        for key, param in self.params.items():
            params[key] = param.value
        return params

@dataclass
class Param():
    key: str
    name: str
    default: any
    type: type
    description: str = ""
    allowed_values: list = None
    nullable: bool = False
    _min: int = None
    _max: int = None
    step: int = None

    def __post_init__(self):
        self.value = self.default
        if self.type == str:
            self.input_type = "text"
        elif self.type == int:
            self.input_type = "number"
        elif self.type == float:
            self.input_type = "number"
        else:
            self.input_type = "text"

    