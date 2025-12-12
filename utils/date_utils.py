"""
date_utils.py
-------------

Utility functions for working with timestamps and datetime objects in UTC. 
Provides helpers for generating formatted timestamps, logging timestamps, 
compact timestamp strings, and converting ISO-formatted strings to datetime objects.

Main Features
-------------
1. Current timestamp
   - now_utc() -> datetime
       Return the current UTC datetime object.

2. String formatting
   - timestamp() -> str
       Return the current UTC timestamp as a string in "YYYY-MM-DD HH:MM:SS UTC" format.
   - ts_compact() -> str
       Return the current UTC timestamp in compact "YYYYMMDD_HHMMSS" format.

3. Logging
   - log_timestamp() -> str
       Return the current UTC timestamp formatted for log messages: "[YYYY-MM-DD HH:MM:SS]".

4. ISO time conversion
   - parse_iso(ts: str) -> datetime
       Convert an ISO-formatted string (optionally ending with "Z") to a UTC datetime object.

Dependencies
------------
- datetime, timezone

Usage Example
-------------
from utils.date_utils import timestamp, ts_compact, log_timestamp, parse_iso

print(timestamp())          # e.g., 2025-12-11 23:45:12 UTC
print(ts_compact())         # e.g., 20251211_234512
print(log_timestamp())      # e.g., [2025-12-11 23:45:12]
dt = parse_iso("2025-12-11T23:45:12Z")
"""


from datetime import datetime, timezone

# -------------------------
# Current timestamp
# -------------------------
def now_utc() -> datetime:
    return datetime.now(timezone.utc)


# -------------------------
# String format
# -------------------------
def timestamp() -> str:
    return now_utc().strftime("%Y-%m-%d %H:%M:%S UTC")


def ts_compact() -> str:
    return now_utc().strftime("%Y%m%d_%H%M%S")


# -------------------------
# For logs
# -------------------------
def log_timestamp() -> str:
    return now_utc().strftime("[%Y-%m-%d %H:%M:%S]")


# -------------------------
# Convert ISO time
# -------------------------
def parse_iso(ts: str) -> datetime:
    return datetime.fromisoformat(ts.replace("Z", "+00:00"))
