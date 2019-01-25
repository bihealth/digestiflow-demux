"""Model with exception classes for the project."""


class DemuxAppException(Exception):
    """Base class for exceptions"""


class ApiProblemException(DemuxAppException):
    """Raised when there is a problem with any API call"""


class MissingConfiguration(DemuxAppException):
    """Raised on missing configuration"""


class InvalidConfiguration(DemuxAppException):
    """Raised on invalid configuration"""


class MissingOutputFile(DemuxAppException):
    """Raised on ``--only-post-message`` and missing output file"""
