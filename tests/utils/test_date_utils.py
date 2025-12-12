# tests/utils/test_date_utils.py

import re
from datetime import datetime, timezone
from utils import date_utils


def test_now_utc_returns_timezone_aware_datetime():
    ts = date_utils.now_utc()
    assert isinstance(ts, datetime)
    assert ts.tzinfo is not None
    assert ts.tzinfo == timezone.utc


def test_timestamp_format():
    ts = date_utils.timestamp()
    # Example: "2025-12-11 20:15:30 UTC"
    pattern = r"^\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2} UTC$"
    assert re.match(pattern, ts)


def test_ts_compact_format():
    ts = date_utils.ts_compact()
    # Example: "20251211_201530"
    pattern = r"^\d{8}_\d{6}$"
    assert re.match(pattern, ts)


def test_log_timestamp_format():
    ts = date_utils.log_timestamp()
    # Example: "[2025-12-11 20:15:30]"
    pattern = r"^\[\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2}\]$"
    assert re.match(pattern, ts)


def test_parse_iso_basic():
    iso = "2025-12-11T20:15:30Z"
    dt = date_utils.parse_iso(iso)
    assert isinstance(dt, datetime)
    assert dt.year == 2025
    assert dt.month == 12
    assert dt.day == 11
    assert dt.hour == 20
    assert dt.minute == 15
    assert dt.second == 30
    assert dt.tzinfo == timezone.utc


def test_parse_iso_with_offset():
    iso = "2025-12-11T20:15:30+00:00"
    dt = date_utils.parse_iso(iso)
    assert dt.tzinfo == timezone.utc
    assert dt.isoformat() == "2025-12-11T20:15:30+00:00"

