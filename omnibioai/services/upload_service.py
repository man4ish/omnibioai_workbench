"""
Module: upload_service
Author: Manish Kumar
Version: 1.0
Date: 2025-12-12

Description:
    Provides the UploadService class for handling file uploads in OmniBioAI.
    Supports validation of file extensions and saving uploaded files to a designated directory.
    Ensures only allowed file types are accepted for security and consistency.

Usage:
    from omnibioai.services.upload_service import UploadService

    service = UploadService(upload_dir="data/uploads")

    # Validate and save a file
    file_content = b"Sample content of the file"
    filename = "example.csv"
    service.validate_file(filename, file_content)
    path = service.save_file(filename, file_content)
    print(f"File saved at: {path}")

Classes:
    - UploadService:
        Service for validating and saving uploaded files.
        
        Methods:
            * __init__(upload_dir="data/uploads"):
                Initializes the service and ensures the upload directory exists.
            * validate_file(filename: str, content: bytes) -> bool:
                Validates that the file extension is allowed. Raises ValueError if unsupported.
            * save_file(filename: str, content: bytes) -> str:
                Saves the validated file to the upload directory and returns the file path.

Constants:
    - ALLOWED_EXTENSIONS: List of file extensions allowed for upload.

Dependencies:
    - os: For directory and file path management.
    - omnibioai.services.logger_service.logger: For logging file operations.
"""


import os
from omnibioai.services.logger_service import logger  # centralized logger

ALLOWED_EXTENSIONS = ['vcf', 'h5', 'csv', 'json', 'pdf']

class UploadService:
    def __init__(self, upload_dir="data/uploads"):
        self.upload_dir = upload_dir
        os.makedirs(upload_dir, exist_ok=True)

    def validate_file(self, filename, content):
        ext = filename.split('.')[-1].lower()
        if ext not in ALLOWED_EXTENSIONS:
            raise ValueError(f"Unsupported file type: {ext}")
        return True

    def save_file(self, filename, content):
        self.validate_file(filename, content)
        path = os.path.join(self.upload_dir, filename)
        with open(path, "wb") as f:
            f.write(content)
        logger.info(f"Uploaded file saved: {path}")
        return path

