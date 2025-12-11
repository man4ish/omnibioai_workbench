#!/bin/bash
# wait-for-service.sh
# Usage: ./wait-for-service.sh <host> <port> [timeout_seconds] -- <command>

set -e

# -----------------------------
# Parse arguments
# -----------------------------
host="$1"
port="$2"
timeout="${3:-60}"  # Default timeout 60 seconds
shift 3
cmd="$@"

if [ -z "$host" ] || [ -z "$port" ] || [ -z "$cmd" ]; then
  echo "Usage: $0 <host> <port> [timeout_seconds] -- <command>"
  exit 1
fi

echo "Waiting for service at $host:$port (timeout: $timeout seconds)..."

# -----------------------------
# Wait loop
# -----------------------------
elapsed=0
while ! nc -z "$host" "$port"; do
  sleep 1
  elapsed=$((elapsed+1))
  echo -n "."
  if [ "$elapsed" -ge "$timeout" ]; then
    echo ""
    echo "Error: Timeout reached while waiting for $host:$port"
    exit 1
  fi
done

echo ""
echo "Service at $host:$port is up! Executing command..."
exec $cmd
