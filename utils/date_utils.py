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
