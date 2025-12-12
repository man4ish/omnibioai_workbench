# test_upload_service.py

import pandas as pd
from omnibioai.services.upload_service import UploadService
from omnibioai.services.logger_service import logger

# ----------------------------
# 1. Upload Example File
# ----------------------------
upload_service = UploadService()
file_path = upload_service.save_file(
    "example.vcf", 
    b"##fileformat=VCF\n#CHROM\tPOS\tID\nchr1\t123\t.\n"
)
logger.info(f"Uploaded file path: {file_path}")

