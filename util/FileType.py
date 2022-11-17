from enum import Enum

class CSV():
    value = 0
    desc = "Comma-separated values"
    extensions = [".csv",".tsv"]

class H5():
    value = 1
    desc = "Hierarchical Data Format"
    extensions = [".h5",".hdf5"]

class RDS():
    value = 2
    desc = "R Data Serialization"
    extensions = [".rds"]

sc_file_types = [CSV, H5, RDS]
sc_file_extensions = [ext for ft in sc_file_types for ext in ft.extensions]
