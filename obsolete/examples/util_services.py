# How Services Will Use These
# upload_service.py
from omnibioai.utils.validators import validate_file
from omnibioai.utils.file_utils import sha256sum

def process_upload(path):
    validate_file(path)
    return sha256sum(path)

# reporting_service.py
from omnibioai.utils.date_utils import timestamp

# rag_service.py
from omnibioai.utils.file_utils import safe_filename