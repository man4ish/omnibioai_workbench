"""
file_utils.py
-------------

Utility functions for common file operations, including hashing, file size calculation, 
safe filename creation, temporary working directories, file copying, and directory creation.

Main Features
-------------
1. Hashing
   - sha256sum(path: str) -> str
       Compute SHA256 checksum of a file.

2. File size
   - get_file_size_mb(path: str) -> float
       Return the size of a file in megabytes.

3. Safe filenames
   - safe_filename(name: str) -> str
       Generate a safe filename by removing invalid characters.

4. Temporary working directory
   - temp_workdir()
       Context manager to create and clean up a temporary working directory.

5. Copy files
   - copy_file(src: str, dst: str)
       Copy a file to a destination, creating parent directories if needed.

6. Directory management
   - ensure_dir(path: str)
       Create a directory if it does not exist.

Dependencies
------------
- hashlib, os, shutil, tempfile, pathlib

Usage Example
-------------
from utils.file_utils import sha256sum, get_file_size_mb, safe_filename, temp_workdir, copy_file, ensure_dir

# Compute SHA256 of a file
checksum = sha256sum("data/file.txt")

# Get file size in MB
size = get_file_size_mb("data/file.txt")

# Use a temporary working directory
with temp_workdir() as d:
    print(f"Temporary dir: {d}")
"""

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
