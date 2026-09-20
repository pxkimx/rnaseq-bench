#!/bin/bash
# Run the server from source with auto-reload (for development with Claude Code).
# Uses the same package environment the app installed: ~/Library/Application Support/RNAseqBench/venv
cd "$(dirname "$0")"
VENV="$HOME/Library/Application Support/RNAseqBench/venv"
[ -x "$VENV/bin/python" ] || VENV=".venv"
[ -x "$VENV/bin/python" ] || { echo "no environment yet — launch the app once (or: python3 -m venv .venv && .venv/bin/pip install -r requirements.txt)"; exit 1; }
export TL_HOME="${TL_HOME:-$HOME/Library/Application Support/RNAseqBench}"
PORT="${PORT:-8765}"
if curl -s -m 1 "http://localhost:$PORT/api/settings" >/dev/null 2>&1; then
  echo "stopping the server already on :$PORT"; curl -s -X POST "http://localhost:$PORT/api/quit" >/dev/null; sleep 1; lsof -ti tcp:"$PORT" | xargs kill 2>/dev/null
fi
echo "http://localhost:$PORT  (auto-reloads when server/ or web/ change; Ctrl+C to stop)"
exec "$VENV/bin/python" -m uvicorn server.app:app --port "$PORT" --reload --reload-dir server --reload-dir web
