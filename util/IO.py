import scanpy as sc
import util.FileFormat as FileFormat

def load_file(path, file_format, file_type, logger):
    if file_format == FileFormat.CSV.value:
        return sc.read_csv(path, delimiter=",")
    elif file_format == FileFormat.TSV.value:
        return sc.read_csv(path, delimiter="\t")
    elif file_format == FileFormat.H5.value:
        return sc.read_h5ad(path)
    elif file_format == FileFormat.RDS.value:
        return sc.read_rds(path)



    