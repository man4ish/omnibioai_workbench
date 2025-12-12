import os
import tempfile
import gzip
import pandas as pd
import pytest
from pathlib import Path
from utils import file_data_utils as fdu


def test_ensure_dir(tmp_path):
    path = tmp_path / "nested/dir"
    out = fdu.ensure_dir(str(path))
    assert os.path.isdir(out)


def test_is_file_and_is_dir(tmp_path):
    file_path = tmp_path / "test.txt"
    file_path.write_text("hello")
    dir_path = tmp_path / "subdir"
    dir_path.mkdir()
    
    assert fdu.is_file(str(file_path))
    assert not fdu.is_file(str(dir_path))
    assert fdu.is_dir(str(dir_path))
    assert not fdu.is_dir(str(file_path))


def test_temp_workdir():
    with fdu.temp_workdir() as d:
        assert os.path.isdir(d)
        tmp_file = Path(d) / "test.txt"
        tmp_file.write_text("data")
        assert tmp_file.exists()
    # directory should be removed after context
    assert not os.path.exists(d)


def test_safe_filename():
    name = "a<>b?c|d.txt"
    safe = fdu.safe_filename(name)
    assert safe == "abcd.txt"


def test_get_file_size_mb(tmp_path):
    file_path = tmp_path / "file.txt"
    file_path.write_bytes(b"0" * 1024 * 1024)  # 1 MB
    size = fdu.get_file_size_mb(str(file_path))
    assert 0.99 <= size <= 1.01


def test_sha256sum(tmp_path):
    file_path = tmp_path / "file.txt"
    file_path.write_bytes(b"hello")
    checksum = fdu.sha256sum(str(file_path))
    assert len(checksum) == 64


def test_copy_file(tmp_path):
    src = tmp_path / "src.txt"
    src.write_text("hello")
    dst = tmp_path / "subdir/dst.txt"
    out = fdu.copy_file(str(src), str(dst))
    assert os.path.isfile(out)
    assert (tmp_path / "subdir/dst.txt").read_text() == "hello"


def test_read_and_save_table(tmp_path):
    df = pd.DataFrame({"a": [1,2], "b": [3,4]})
    path_csv = tmp_path / "data.csv"
    fdu.save_table(df, str(path_csv))
    df_read = fdu.read_table(str(path_csv))
    pd.testing.assert_frame_equal(df, df_read)
    
    path_excel = tmp_path / "data.xlsx"
    fdu.save_table(df, str(path_excel), excel=True)
    df_read2 = fdu.read_table(str(path_excel), excel=True)
    pd.testing.assert_frame_equal(df, df_read2)


def test_open_file(tmp_path):
    file_path = tmp_path / "test.txt"
    file_path.write_text("hello")
    with fdu.open_file(str(file_path)) as f:
        assert f.read().strip() == "hello"

    gz_path = tmp_path / "test.gz"
    with gzip.open(gz_path, "wt") as f:
        f.write("world")
    with fdu.open_file(str(gz_path)) as f:
        assert f.read().strip() == "world"

