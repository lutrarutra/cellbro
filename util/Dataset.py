import util.FileFormat as FileFormat, util.FileType as FileType, util.IO as IO


class Dataset():
    def __init__(self, path, file_type, logger):
        self.path = path
        self.file_type = file_type

        if path.endswith(".csv"):
            self.file_format = FileFormat.CSV.value
        elif path.endswith(".tsv"):
            self.file_format = FileFormat.TSV.value
        elif path.endswith(".h5") or path.endswith(".hdf5"):
            self.file_format = FileFormat.H5.value
        elif path.endswith(".rds"):
            self.file_format = FileFormat.RDS.value
        else:
            logger.error(f"File format {'.' + path.split('.')[-1]} not supported...")
            return None

        self.logger = logger
        self.adata = IO.load_file(path, self.file_format, self.file_type, self.logger)
        