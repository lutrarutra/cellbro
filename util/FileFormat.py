class CSV():
    value = 0
    desc = "Comma-separated values"
    extensions = [".csv"]

class TSV():
    value = 1
    desc = "Tab-separated values"
    extensions = [".tsv"]

class H5():
    value = 2
    desc = "Hierarchical Data Format"
    extensions = [".h5",".hdf5"]

class RDS():
    value = 3
    desc = "R Data Serialization"
    extensions = [".rds"]

sc_file_formats = [CSV, TSV, H5, RDS]
sc_file_extensions = [ext for ft in sc_file_formats for ext in ft.extensions]
