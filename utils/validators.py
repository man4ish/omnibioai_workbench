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
