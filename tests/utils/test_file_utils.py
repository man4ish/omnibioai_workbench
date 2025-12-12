# tests/utils/test_file_utils.py

import os
import hashlib
import tempfile
from pathlib import Path
from contextlib import contextmanager

import pytest

from utils.file_utils import (
    sha256sum,
    get_file_size_mb,
    safe_filename,
    temp_workdir,
    copy_file,
    ensure_dir,
)


# -------------------------------------------------------------------
# Helper to convert the generator into a real context manager
# -------------------------------------------------------------------
@contextmanager
def tempdir_cm():
    gen = temp_workdir()
    d = next(gen)  # get directory from generator
    try:
        yield d
    finally:
        try:
            next(gen)
        except StopIteration:
            pass


# -------------------------------------------------------------------
# Test SHA256 hashing
# -------------------------------------------------------------------
def test_sha256sum(tmp_path):
    test_file = tmp_path / "sample.txt"
    test_file.write_text("hello world")

    expected = hashlib.sha256(b"hello world").hexdigest()
    result = sha256sum(str(test_file))

    assert result == expected


# -------------------------------------------------------------------
# Test file size in MB
# -------------------------------------------------------------------
def test_get_file_size_mb(tmp_path):
    test_file = tmp_path / "small.bin"
    test_file.write_bytes(b"a" * 1024 * 1024)  # 1 MB

    size = get_file_size_mb(str(test_file))
    assert 0.99 < size < 1.01  # allow slight FS rounding


# -------------------------------------------------------------------
# Test safe filename
# -------------------------------------------------------------------
def test_safe_filename():
    name = "a/b*c?d e:f|g.txt"
    cleaned = safe_filename(name)
    assert cleaned == "abcdefg.txt"  # removes illegal characters


# -------------------------------------------------------------------
# Test temp_workdir context manager
# -------------------------------------------------------------------
def test_temp_workdir_creates_and_deletes():
    with tempdir_cm() as d:
        d_path = Path(d)
        assert d_path.exists()
        assert d_path.is_dir()

        # create a file inside it
        file_inside = d_path / "test.txt"
        file_inside.write_text("test ok")
        assert file_inside.exists()

    # After exiting, directory must be removed
    assert not Path(d).exists()


# -------------------------------------------------------------------
# Test copy file
# -------------------------------------------------------------------
def test_copy_file(tmp_path):
    src = tmp_path / "src.txt"
    src.write_text("copy me")

    dst = tmp_path / "subdir" / "dst.txt"

    final = copy_file(str(src), str(dst))

    assert final == str(dst)
    assert dst.exists()
    assert dst.read_text() == "copy me"


# -------------------------------------------------------------------
# Test ensure_dir
# -------------------------------------------------------------------
def test_ensure_dir(tmp_path):
    new_dir = tmp_path / "mydir"
    ensure_dir(str(new_dir))

    assert new_dir.exists()
    assert new_dir.is_dir()

