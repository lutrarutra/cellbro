from abc import ABC, abstractmethod


class FileFormat(ABC):
    @staticmethod
    @abstractmethod
    def ext():
        ...

    @staticmethod
    @abstractmethod
    def name():
        ...

    @staticmethod
    @abstractmethod
    def desc():
        ...

class CSV(FileFormat):
    @staticmethod
    def ext():
        return ".csv"

    @staticmethod
    def name():
        return "CSV"

    @staticmethod
    def desc():
        return "Comma Separated Values"

class TSV(FileFormat):
    @staticmethod
    def ext():
        return ".tsv"

    @staticmethod
    def name():
        return "TSV"

    @staticmethod
    def desc():
        return "Tab Separated Values"

class Pickle(FileFormat):
    @staticmethod
    def ext():
        return ".pkl"

    @staticmethod
    def name():
        return "Pickle"

    @staticmethod
    def desc():
        return "Binary"


class H5AD(FileFormat):
    @staticmethod
    def ext():
        return ".h5ad"

    @staticmethod
    def name():
        return "H5 AnnData"

    @staticmethod
    def desc():
        return "Hierarchical Data Format"




