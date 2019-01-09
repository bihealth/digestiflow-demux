"""Model with exception classes for the project."""


class DemuxAppException(Exception):
    """Base class for exceptions"""


class MissingConfiguration(DemuxAppException):
    """Raised on missing configuration"""


class InvalidConfiguration(DemuxAppException):
    """Raised on invalid configuration"""
