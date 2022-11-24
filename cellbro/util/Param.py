
class ParamsDict():
    def __init__(self, params_list):
        self.params = {}
        for param in params_list:
            self.params[param.key] = param

    def update(self, params_dict):
        for i, key in enumerate(list(self.keys())):
            self.params[key].value = params_dict[key]

        return self

    def __getitem__(self, key):
        return self.params[key]

    def add(self, param):
        self.params[param.key] = param

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

class Param():
    def __init__(
        self, key, name, default, type,
        description, nullable=False,
        _min=None, _max=None, step=None,
        ):
        self.key = key
        self.name = name
        self.default = default
        self.value = default
        self.type = type
        if self.type == str:
            self.input_type = "text"
        elif self.type == int:
            self.input_type = "number"
        elif self.type == float:
            self.input_type = "number"
        else:
            self.input_type = "text"
        self.description = description
        self.nullable = nullable
        self.min = _min
        self.max = _max
        self.step = step

    