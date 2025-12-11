import hashlib
import os
import shutil
import tempfile
from pathlib import Path


# -------------------------
# Hashing (SHA256)
# -------------------------
def sha256sum(path: str) -> str:
    h = hashlib.sha256()

    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)

    return h.hexdigest()


# -------------------------
# File size
# -------------------------
def get_file_size_mb(path: str) -> float:
    return os.path.getsize(path) / (1024 * 1024)


# -------------------------
# Create safe filename
# -------------------------
def safe_filename(name: str) -> str:
    return "".join(c for c in name if c.isalnum() or c in ("-", "_", ".")).strip()


# -------------------------
# Temporary working directory
# -------------------------
def temp_workdir():
    """
    Usage:
        with temp_workdir() as d:
            # write files inside d
    """
    d = tempfile.mkdtemp()
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)


# -------------------------
# Copy with overwrite
# -------------------------
def copy_file(src: str, dst: str):
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)
    return dst


# -------------------------
# Ensure directory exists
# -------------------------
def ensure_dir(path: str):
    Path(path).mkdir(parents=True, exist_ok=True)
    return path
