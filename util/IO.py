import scanpy as sc
import util.FileFormat as FileFormat

def load_file(path, file_format, file_type, logger):
    if file_format.value == FileFormat.CSV.value:
        return sc.read_csv(path, delimiter=",")
    elif file_format.value == FileFormat.TSV.value:
        return sc.read_csv(path, delimiter="\t")
    elif file_format.value == FileFormat.H5.value:
        return sc.read_h5ad(path)
    elif file_format.value == FileFormat.RDS.value:
        return sc.read_rds(path)



    