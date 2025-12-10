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

