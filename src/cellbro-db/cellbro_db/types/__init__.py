from enum import IntEnum

class UserType(IntEnum):
    DEACTIVATED = 0
    REGULAR = 10
    ADMIN = 100

class PaperProcessingStatus(IntEnum):
    DRAFT = 10
    PROCESSING = 20
    PROCESSED = 100
    MISSING_METADATA = 200
    INVALID_PDF_URL = 201
    ERROR = 300


class SourceType(IntEnum):
    UNKNOWN = 0
    REPOSITORY = 1
    JOURNAL = 2

class PaperType(IntEnum):
    UNKNOWN = 0
    PREPRINT = 1
    ARTICLE = 2
    DATASET = 3

class OADomain(IntEnum):
    UNKNOWN = 0
    LIFE_SCIENCES = 1
    SOCIAL_SCIENCES = 2
    PHYSICAL_SCIENCES = 3
    HEALTH_SCIENCES = 4

class AuthorPosition(IntEnum):
    FIRST = 1
    MIDDLE = 2
    LAST = 3
    UNSPECIFIED = 100

