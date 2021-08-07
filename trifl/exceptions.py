class TriflError(Exception):
    """Base class for exceptions in this project."""

class DtaFileNotFound(TriflError):
    """Raised when a DTASelect file is not found at the expected path."""

class CimageFlatFileNotFound(TriflError):
    """Raised when a cimage flatfile (output_to_excel.txt) is not found at the expected path."""

class FilterNotFoundException(TriflError):
    """Raised when filter is not found."""
