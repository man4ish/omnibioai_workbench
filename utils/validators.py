"""
validators.py
-------------

Utility functions for validating different file types including VCF, CSV, HDF5, and JSON. 
Raises `ValidationError` for invalid files or unsupported formats.

Main Features
-------------
1. Custom exception
   - ValidationError
       Exception raised when a file fails validation.

2. VCF validation
   - validate_vcf(path: str) -> bool
       Validate a VCF file (.vcf or .vcf.gz) by checking existence and mandatory header (#CHROM).

3. CSV validation
   - validate_csv(path: str, min_columns: int = 1) -> bool
       Validate a CSV file, ensuring it exists and has at least `min_columns` in the header.

4. HDF5 validation
   - validate_h5(path: str) -> bool
       Validate an HDF5 (.h5 or .hdf5) file by opening the root group.

5. JSON validation
   - validate_json(path: str) -> bool
       Validate a JSON file by attempting to load it.

6. File type router
   - validate_file(path: str) -> bool
       Dispatches validation based on file extension. Supports VCF, CSV, HDF5, and JSON.

Dependencies
------------
- os, gzip, h5py, csv, json

Usage Example
-------------
from utils.validators import validate_file, ValidationError

try:
    validate_file("data/sample.vcf.gz")
except ValidationError as e:
    print(f"Validation failed: {e}")
"""

import os
import gzip
import h5py
import csv
import json

class ValidationError(Exception):
    pass


# -------------------------
# VCF Validation
# -------------------------
def validate_vcf(path: str) -> bool:
    if not os.path.exists(path):
        raise ValidationError("VCF file does not exist.")

    opener = gzip.open if path.endswith(".gz") else open

    try:
        with opener(path, "rt") as f:
            for line in f:
                if line.startswith("#CHROM"):
                    return True
    except Exception as e:
        raise ValidationError(f"Invalid VCF file: {e}")

    raise ValidationError("VCF missing mandatory header (#CHROM).")


# -------------------------
# CSV Validation
# -------------------------
def validate_csv(path: str, min_columns: int = 1) -> bool:
    if not os.path.exists(path):
        raise ValidationError("CSV file does not exist.")

    try:
        with open(path, newline="") as f:
            reader = csv.reader(f)
            header = next(reader, None)

            if not header or len(header) < min_columns:
                raise ValidationError("CSV header invalid or too few columns.")

        return True

    except Exception as e:
        raise ValidationError(f"Invalid CSV file: {e}")


# -------------------------
# HDF5 Validation
# -------------------------
def validate_h5(path: str) -> bool:
    if not os.path.exists(path):
        raise ValidationError("H5 file does not exist.")

    try:
        with h5py.File(path, "r") as f:
            # Check if root can be opened
            _ = list(f.keys())
        return True

    except Exception as e:
        raise ValidationError(f"Invalid HDF5 file: {e}")


# -------------------------
# JSON Validation
# -------------------------
def validate_json(path: str) -> bool:
    if not os.path.exists(path):
        raise ValidationError("JSON file does not exist.")

    try:
        json.load(open(path))
        return True
    except Exception as e:
        raise ValidationError(f"Invalid JSON file: {e}")


# -------------------------
# File type router
# -------------------------
def validate_file(path: str) -> bool:
    """
    Dispatch validator based on file type.
    """
    if path.endswith(".vcf") or path.endswith(".vcf.gz"):
        return validate_vcf(path)

    if path.endswith(".csv"):
        return validate_csv(path)

    if path.endswith(".h5") or path.endswith(".hdf5"):
        return validate_h5(path)

    if path.endswith(".json"):
        return validate_json(path)

    raise ValidationError(f"Unsupported file type: {path}")
