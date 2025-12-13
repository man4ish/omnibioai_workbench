"""
exceptions.py

Custom exception classes for OmnibioAI services.

This module defines specialized exceptions used throughout
the OmnibioAI codebase to improve error handling and clarity.

Exceptions
----------
DatasetNotFoundError
    Raised when a requested dataset is not found in the local registry or cache.

AnnotationError
    Raised when an error occurs during variant annotation or data processing.

MLModelNotFoundError
    Raised when a requested machine learning model is not available or cannot be located.
"""

class DatasetNotFoundError(Exception):
    pass

class AnnotationError(Exception):
    pass

class MLModelNotFoundError(Exception):
    pass
