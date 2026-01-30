from enum import IntEnum, StrEnum

class TaskStatus(StrEnum):
    RUNNING = "running"
    COMPLETED = "completed"
    FAILED = "failed"

class UserType(IntEnum):
    DEACTIVATED = 0
    REGULAR = 10
    ADMIN = 100

class ChecklistStep(IntEnum):
    LOAD = 1
    AGGREATE = 2
    QC = 3
    NORMALIZE = 4
    DIM_REDUCTION = 5
    NEIHBOR_GRAPH_CONSTRUCTION = 6
    CLUSTERING = 7
    PROJECTION = 8
    DEA = 9
    GSEA = 10


class VariableType(IntEnum):
    STRING = 1
    CATEGORICAL = 2
    INTEGER = 3
    FLOAT = 4
    BOOLEAN = 5
    LIST = 6
    ARRAY = 7
    MATRIX = 8

class AnnDataLayerType(StrEnum):
    OBS = "obs"
    VAR = "var"
    UNS = "uns"
    LAYER = "layer"
    OBSM = "obsm"
    OBSP = "obsp"
    VARM = "varm"
    VARP = "varp"
    X = "X"


