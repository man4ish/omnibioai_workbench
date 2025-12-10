#!/bin/bash

# wait-for-ollama.sh

set -e

host="$1"
port="$2"
shift 2
cmd="$@"

echo "Waiting for Ollama server at $host:$port..."

while ! nc -z "$host" "$port"; do
  sleep 1
done

echo "Ollama server is up - executing command"
exec $cmd

