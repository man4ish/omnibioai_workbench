import os
import gzip
import shutil
import tempfile
import hashlib
import pandas as pd
from pathlib import Path
import requests
from contextlib import contextmanager

# -------------------------
# File existence & path helpers
# -------------------------
def ensure_dir(path: str) -> str:
    """Create directory if it does not exist."""
    Path(path).mkdir(parents=True, exist_ok=True)
    return path

def is_file(path: str) -> bool:
    """Check if a path exists and is a file."""
    return os.path.isfile(path)

def is_dir(path: str) -> bool:
    """Check if a path exists and is a directory."""
    return os.path.isdir(path)


# -------------------------
# Temporary working directory
# -------------------------
@contextmanager
def temp_workdir():
    """Context manager for a temporary working directory."""
    d = Path(tempfile.mkdtemp())
    try:
        yield d
    finally:
        shutil.rmtree(d, ignore_errors=True)


# -------------------------
# Safe filename
# -------------------------
def safe_filename(name: str) -> str:
    """Make a filename safe by removing invalid characters."""
    return "".join(c for c in name if c.isalnum() or c in ("-", "_", ".")).strip()


# -------------------------
# File size
# -------------------------
def get_file_size_mb(path: str) -> float:
    """Return file size in megabytes."""
    return os.path.getsize(path) / (1024 * 1024)


# -------------------------
# Hashing (SHA256)
# -------------------------
def sha256sum(path: str) -> str:
    """Compute SHA256 checksum of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


# -------------------------
# Copy with overwrite
# -------------------------
def copy_file(src: str, dst: str) -> str:
    """Copy file, creating parent directories if needed."""
    Path(dst).parent.mkdir(parents=True, exist_ok=True)
    shutil.copy(src, dst)
    return dst


# -------------------------
# Read CSV / TSV / Excel
# -------------------------
def read_table(path: str, sep="\t", excel=False) -> pd.DataFrame:
    """Read CSV/TSV or Excel file into pandas DataFrame."""
    if not is_file(path):
        raise FileNotFoundError(f"{path} not found")
    
    if excel:
        return pd.read_excel(path)
    else:
        return pd.read_csv(path, sep=sep)


# -------------------------
# Write CSV / TSV / Excel
# -------------------------
def save_table(df: pd.DataFrame, path: str, sep="\t", excel=False) -> str:
    """Save pandas DataFrame to CSV/TSV or Excel."""
    ensure_dir(os.path.dirname(path))
    if excel:
        df.to_excel(path, index=False)
    else:
        df.to_csv(path, sep=sep, index=False)
    return path


# -------------------------
# Open gzip or plain text
# -------------------------
def open_file(path: str, mode="rt"):
    """Open plain text or gzipped file."""
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)

