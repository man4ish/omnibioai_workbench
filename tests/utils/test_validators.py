import os
import gzip
import json
import h5py
import pytest
from utils.validators import (
    validate_vcf,
    validate_csv,
    validate_h5,
    validate_json,
    validate_file,
    ValidationError,
)

# -------------------------
# Helpers
# -------------------------

def write(path, text):
    with open(path, "w") as f:
        f.write(text)


# -------------------------
# VCF TESTS
# -------------------------

def test_validate_vcf_valid(tmp_path):
    p = tmp_path / "good.vcf"
    write(str(p), "##fileformat=VCFv4.2\n#CHROM\tPOS\n1\t100\n")
    assert validate_vcf(str(p)) is True


def test_validate_vcf_valid_gz(tmp_path):
    p = tmp_path / "good.vcf.gz"
    with gzip.open(p, "wt") as f:
        f.write("##meta\n#CHROM\tPOS\n1\t100\n")
    assert validate_vcf(str(p)) is True


def test_validate_vcf_missing_header(tmp_path):
    p = tmp_path / "bad.vcf"
    write(str(p), "##meta\n1\t100\n")
    with pytest.raises(ValidationError):
        validate_vcf(str(p))


def test_validate_vcf_missing_file(tmp_path):
    p = tmp_path / "missing.vcf"
    with pytest.raises(ValidationError):
        validate_vcf(str(p))


# -------------------------
# CSV TESTS
# -------------------------

def test_validate_csv_valid(tmp_path):
    p = tmp_path / "good.csv"
    write(str(p), "col1,col2\n1,2\n")
    assert validate_csv(str(p))


def test_validate_csv_invalid_header(tmp_path):
    p = tmp_path / "bad.csv"
    write(str(p), "\n1,2\n")
    with pytest.raises(ValidationError):
        validate_csv(str(p))


def test_validate_csv_missing(tmp_path):
    p = tmp_path / "missing.csv"
    with pytest.raises(ValidationError):
        validate_csv(str(p))


# -------------------------
# HDF5 TESTS
# -------------------------

def test_validate_h5_valid(tmp_path):
    p = tmp_path / "good.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("x", data=[1, 2, 3])

    assert validate_h5(str(p))


def test_validate_h5_invalid(tmp_path):
    p = tmp_path / "bad.h5"
    write(str(p), "not an hdf5 file at all")
    with pytest.raises(ValidationError):
        validate_h5(str(p))


def test_validate_h5_missing(tmp_path):
    p = tmp_path / "missing.h5"
    with pytest.raises(ValidationError):
        validate_h5(str(p))


# -------------------------
# JSON TESTS
# -------------------------

def test_validate_json_valid(tmp_path):
    p = tmp_path / "good.json"
    with open(p, "w") as f:
        json.dump({"a": 1}, f)
    assert validate_json(str(p))


def test_validate_json_invalid(tmp_path):
    p = tmp_path / "bad.json"
    write(str(p), "{not json}")
    with pytest.raises(ValidationError):
        validate_json(str(p))


def test_validate_json_missing(tmp_path):
    p = tmp_path / "missing.json"
    with pytest.raises(ValidationError):
        validate_json(str(p))


# -------------------------
# Dispatcher TESTS
# -------------------------

def test_validate_file_dispatch_vcf(tmp_path):
    p = tmp_path / "file.vcf"
    write(str(p), "##fileformat=VCF\n#CHROM\tPOS\n1\t100")
    assert validate_file(str(p))


def test_validate_file_dispatch_csv(tmp_path):
    p = tmp_path / "file.csv"
    write(str(p), "a,b\n1,2")
    assert validate_file(str(p))


def test_validate_file_dispatch_h5(tmp_path):
    p = tmp_path / "file.h5"
    with h5py.File(p, "w") as f:
        f.create_dataset("d", data=[1])
    assert validate_file(str(p))


def test_validate_file_dispatch_json(tmp_path):
    p = tmp_path / "file.json"
    with open(p, "w") as f:
        json.dump({"x": 1}, f)
    assert validate_file(str(p))


def test_validate_file_unsupported_extension(tmp_path):
    p = tmp_path / "file.txt"
    write(str(p), "hello")
    with pytest.raises(ValidationError):
        validate_file(str(p))

