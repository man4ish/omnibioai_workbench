"""
download_upload.py
------------------

Utility functions for downloading and uploading files over HTTP. 
Supports retries, automatic directory creation, and optional extraction of compressed files.

Main Features
-------------
1. Download utility
   - download_file(url: str, output_dir: str = "data", overwrite: bool = False, unzip: bool = True, retries: int = 3) -> str
       Download a file from a URL to the specified directory.
       Supports retry attempts, overwriting existing files, and automatic extraction for .gz and .zip files.

2. Upload utility
   - upload_file(file_path: str, url: str, retries: int = 3) -> bool
       Upload a file to a specified URL using HTTP POST.
       Supports retry attempts and raises FileNotFoundError if the file does not exist.

Dependencies
------------
- os, gzip, shutil, zipfile, pathlib, requests

Usage Example
-------------
from utils.download_upload import download_file, upload_file

# Download a file and automatically unzip if compressed
path = download_file("https://example.com/data.zip")

# Upload a file
success = upload_file("data/file.txt", "https://example.com/upload")
"""

import os
import gzip
import shutil
import zipfile
from pathlib import Path
import requests

# -------------------------
# Download utility
# -------------------------
def download_file(
    url: str,
    output_dir: str = "data",
    overwrite: bool = False,
    unzip: bool = True,
    retries: int = 3
) -> str:
    os.makedirs(output_dir, exist_ok=True)
    filename = url.split("/")[-1]
    filepath = os.path.join(output_dir, filename)

    if os.path.exists(filepath) and not overwrite:
        print(f"{filename} already exists, skipping download.")
    else:
        for attempt in range(retries):
            try:
                print(f"Downloading {filename} (Attempt {attempt+1})...")
                with requests.get(url, stream=True, timeout=30) as r:
                    r.raise_for_status()
                    with open(filepath, 'wb') as f:
                        for chunk in r.iter_content(chunk_size=8192):
                            if chunk:
                                f.write(chunk)
                print(f"Downloaded {filename} to {output_dir}")
                break
            except Exception as e:
                print(f"Download failed (Attempt {attempt+1}): {e}")
                if attempt == retries - 1:
                    raise

    if unzip:
        if filepath.endswith(".gz"):
            extracted_path = os.path.join(output_dir, filename[:-3])
            with gzip.open(filepath, 'rb') as f_in, open(extracted_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
            print(f"Extracted {filename} to {extracted_path}")
            return extracted_path
        elif filepath.endswith(".zip"):
            with zipfile.ZipFile(filepath, 'r') as zip_ref:
                zip_ref.extractall(output_dir)
            print(f"Extracted {filename} into {output_dir}")
            return output_dir

    return filepath

# -------------------------
# Upload utility
# -------------------------
def upload_file(file_path: str, url: str, retries: int = 3) -> bool:
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"{file_path} does not exist.")

    for attempt in range(retries):
        try:
            with open(file_path, 'rb') as f:
                response = requests.post(url, files={"file": f})
                response.raise_for_status()
            print(f"Uploaded {file_path.name} to {url}")
            return True
        except Exception as e:
            print(f"Upload failed (Attempt {attempt+1}): {e}")
            if attempt == retries - 1:
                raise
    return False
