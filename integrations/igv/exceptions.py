class IGVError(Exception):
    """Base class for all IGV-related errors."""
    pass


class IGVConnectionError(IGVError):
    """Raised when IGV socket or IGV application cannot be reached."""
    pass


class IGVCommandError(IGVError):
    """Raised when an IGV command fails or IGV returns an error."""
    pass
