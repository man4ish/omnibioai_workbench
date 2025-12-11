class RAGBioError(Exception):
    """Base class for all RAGBio errors."""
    pass


class RAGBioConnectionError(RAGBioError):
    """Raised when RAGBio services or indexes cannot be reached."""
    pass


class RAGBioQueryError(RAGBioError):
    """Raised when a RAGBio search or embedding fails."""
    pass


class RAGBioPipelineError(RAGBioError):
    """Raised when a pipeline operation fails."""
    pass
