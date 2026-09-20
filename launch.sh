#!/bin/bash
# RNAseq Bench launcher (shared by "RNAseq Bench.app" and "Start RNAseq Bench.command").
#   MODE=gui       -> notifications + a small "running" dialog with a Quit button (used by the .app)
#   MODE=terminal  -> prints everything to the Terminal window and runs the server in the foreground
# Environment lives in ~/Library/Application Support/RNAseqBench (venv, jobs, settings, log).
set -u
# If macOS started us under Rosetta (Intel translation) on an Apple-Silicon Mac, restart natively:
# otherwise pip installs Intel packages, several of which have to be compiled and fail.
if [ "$(sysctl -n sysctl.proc_translated 2>/dev/null)" = "1" ] && [ -z "${RB_NATIVE:-}" ]; then
  export RB_NATIVE=1
  exec arch -arm64 /bin/bash "$0" "$@"
fi
MODE="${MODE:-terminal}"
APP_DIR="$(cd "$(dirname "$0")" && pwd)"
HOME_DIR="$HOME/Library/Application Support/RNAseqBench"
VENV="$HOME_DIR/venv"
LOG="$HOME_DIR/rnaseq-bench.log"
PORT="${PORT:-8765}"
VERSION="$(cat "$APP_DIR/VERSION" 2>/dev/null || echo dev)"
# one-time migration from the old name: analyses and settings come along, the package environment is
# rebuilt (a Python venv cannot be moved — it hard-codes its original path)
OLD_DIR="$HOME/Library/Application Support/TranscriptLens"
if [ ! -d "$HOME_DIR" ] && [ -d "$OLD_DIR" ]; then
  mkdir -p "$HOME_DIR"
  for item in jobs transcript-lens.log; do [ -e "$OLD_DIR/$item" ] && mv "$OLD_DIR/$item" "$HOME_DIR/" 2>/dev/null; done
fi
mkdir -p "$HOME_DIR"
ARCH="$(uname -m)"
export TL_HOME="$HOME_DIR"
export PATH="/opt/homebrew/bin:/usr/local/bin:/Library/Frameworks/Python.framework/Versions/3.13/bin:/Library/Frameworks/Python.framework/Versions/3.12/bin:/Library/Frameworks/Python.framework/Versions/3.11/bin:/Library/Frameworks/Python.framework/Versions/3.10/bin:$PATH"
export PYTHONUNBUFFERED=1

log()  { echo "[$(date '+%H:%M:%S')] $*" >>"$LOG"; [ "$MODE" = terminal ] && echo "$*"; }
note() { log "$*"; [ "$MODE" = gui ] && osascript -e "display notification \"$1\" with title \"RNAseq Bench $VERSION\"" >/dev/null 2>&1; }
fail() {
  log "!! $1"
  if [ "$MODE" = gui ]; then
    osascript -e "display dialog \"$1\n\nLog: $LOG\" with title \"RNAseq Bench $VERSION\" buttons {\"Open log\",\"OK\"} default button \"OK\" with icon stop" 2>/dev/null | grep -q "Open log" && open -a Console "$LOG"
  else
    echo; echo "Log: $LOG"; read -r -p "Press Enter to close."
  fi
  exit 1
}

log "=== RNAseq Bench $VERSION starting ($MODE) — app code: $APP_DIR"

# ---- a server already on the port? reuse it if it is this version, otherwise replace it
RUNNING_VER="$(curl -s -m 2 "http://localhost:$PORT/api/settings" 2>/dev/null | sed -n 's/.*"version": *"\([^"]*\)".*/\1/p')"
if [ "${1:-}" = "--selftest" ]; then
  RUNNING_VER=""          # the self-test never touches the port
elif [ -n "$RUNNING_VER" ]; then
  if [ "$RUNNING_VER" = "$VERSION" ]; then
    log "Already running (v$RUNNING_VER) — opening the browser."
    open "http://localhost:$PORT"
    exit 0
  fi
  log "An older RNAseq Bench (v$RUNNING_VER) is on port $PORT — stopping it."
  lsof -ti tcp:"$PORT" 2>/dev/null | xargs kill 2>/dev/null
  for _ in 1 2 3 4 5 6 7 8 9 10; do curl -s -m 1 "http://localhost:$PORT/api/settings" >/dev/null 2>&1 || break; sleep 1; done
elif lsof -ti tcp:"$PORT" >/dev/null 2>&1; then
  fail "Port $PORT is used by another program. Quit it, or start RNAseq Bench with a different port: PORT=8766 \"$0\""
fi

# ---- find Python 3.10–3.13 (3.14 not yet: harmonypy has no wheel and scanpy import breaks)
PY=""
# Only a Python built for this Mac's chip is acceptable: an Intel build running under Rosetta on an
# Apple-Silicon Mac installs packages that cannot load ("incompatible architecture").
for c in python3.13 python3.12 python3.11 python3.10 python3 \
         /Library/Frameworks/Python.framework/Versions/3.13/bin/python3.13 \
         /Library/Frameworks/Python.framework/Versions/3.12/bin/python3.12 \
         /opt/homebrew/bin/python3.13 /opt/homebrew/bin/python3.12 /opt/homebrew/bin/python3.11 \
         /usr/bin/python3; do
  command -v "$c" >/dev/null 2>&1 || continue
  info=$("$c" -c 'import sys,platform;print(sys.version_info[1], platform.machine())' 2>/dev/null) || continue
  v="${info%% *}"; m="${info##* }"
  if [ "$v" -ge 10 ] && [ "$v" -le 13 ] && [ "$m" = "$ARCH" ]; then PY="$(command -v "$c")"; break; fi
  log "skipping $(command -v "$c") (3.$v, $m; this Mac is $ARCH)"
done
if [ -z "$PY" ]; then
  if [ "$MODE" = gui ]; then
    osascript -e 'display dialog "RNAseq Bench needs Python 3.10–3.13 (3.14 is not supported by the analysis libraries yet).\n\nInstall Python 3.13 from python.org (recommended), or run\n    xcode-select --install\nin Terminal. Then open RNAseq Bench again." with title "RNAseq Bench" buttons {"Open python.org","OK"} default button "OK" with icon caution' 2>/dev/null | grep -q python.org && open "https://www.python.org/downloads/macos/"
    exit 1
  fi
  fail "Python 3.10–3.13 not found (3.14 is not supported by the analysis libraries yet). Install 3.13 from python.org (or: xcode-select --install) and run this again."
fi
log "Python: $PY ($("$PY" --version 2>&1))"

# ---- (re)install packages when the environment is missing or requirements changed
# requirements.txt is split at the "--- optional" marker: the core block must install for the app to
# run at all, the optional block only powers add-ons that degrade to a note, so it is installed
# best-effort and retried on the next launch rather than blocking startup.
CORE_REQ="$HOME_DIR/.core-requirements.txt"
OPT_REQ="$HOME_DIR/.optional-requirements.txt"
mkdir -p "$HOME_DIR"
awk '/^# --- optional/{o=1} !o' "$APP_DIR/requirements.txt" >"$CORE_REQ"
awk '/^# --- optional/{o=1} o'  "$APP_DIR/requirements.txt" >"$OPT_REQ"
REQ_SHA="$(shasum -a 256 "$CORE_REQ" | cut -c1-16)"
OPT_SHA="$(shasum -a 256 "$OPT_REQ" | cut -c1-16)"
# a broken environment (moved folder, wrong chip architecture, deleted base Python) is thrown away and rebuilt
if [ -d "$VENV" ]; then
  vinfo=$("$VENV/bin/python" -c 'import platform,sys;print(platform.machine(), sys.prefix)' 2>/dev/null) || vinfo=""
  # packages installed for the other architecture (from an earlier Rosetta launch) make numpy fail to load
  if [ -n "$vinfo" ] && ! "$VENV/bin/python" -c "import numpy" >/dev/null 2>&1 && [ -d "$VENV/lib" ]; then vinfo=""; fi
  if [ -z "$vinfo" ] || [ "${vinfo%% *}" != "$ARCH" ] || ! "$VENV/bin/pip" --version >/dev/null 2>&1; then
    log "Package environment is unusable (${vinfo:-no python}) — rebuilding it."
    rm -rf "$VENV"
  fi
fi
NEED_INSTALL=0
if [ ! -x "$VENV/bin/python" ]; then NEED_INSTALL=1
elif [ "$(cat "$VENV/.requirements.sha" 2>/dev/null)" != "$REQ_SHA" ]; then NEED_INSTALL=1
elif ! "$VENV/bin/python" -c "import scanpy, pydeseq2, fastapi, reportlab, anthropic, skimage, openpyxl, harmonypy, leidenalg" >/dev/null 2>&1; then NEED_INSTALL=1
fi
NEED_OPTIONAL=0
[ "$(cat "$VENV/.optional.sha" 2>/dev/null)" = "$OPT_SHA" ] || NEED_OPTIONAL=1
if [ "$NEED_INSTALL" = 1 ]; then
  if [ ! -x "$VENV/bin/python" ]; then
    note "First launch: installing analysis packages (3–5 minutes)."
    "$PY" -m venv "$VENV" >>"$LOG" 2>&1 || fail "Could not create a Python environment with $PY."
  else
    note "Updating analysis packages for v$VERSION (1–3 minutes)."
  fi
  if [ "$MODE" = terminal ]; then
    "$VENV/bin/pip" install --upgrade pip 2>&1 | tee -a "$LOG" | grep -v "already satisfied"
    "$VENV/bin/pip" install -r "$CORE_REQ" 2>&1 | tee -a "$LOG" | grep -v "already satisfied"
  else
    "$VENV/bin/pip" install --upgrade pip >>"$LOG" 2>&1
    "$VENV/bin/pip" install -r "$CORE_REQ" >>"$LOG" 2>&1
  fi
  "$VENV/bin/python" -c "import scanpy, pydeseq2, fastapi, reportlab, anthropic, skimage, openpyxl, harmonypy, leidenalg" >>"$LOG" 2>&1 \
    || { rm -f "$VENV/.requirements.sha"; fail "Installing packages failed. Common causes: no internet, or a Python build without SSL. Details are in the log."; }
  echo "$REQ_SHA" >"$VENV/.requirements.sha"
  log "Packages installed."
fi

# ---- expert add-ons: best-effort, never fatal
if [ "$NEED_OPTIONAL" = 1 ]; then
  note "Installing packages for the expert analyses (GSEA, pathway activity, cell-cell communication)."
  if [ "$MODE" = terminal ]; then
    "$VENV/bin/pip" install -r "$OPT_REQ" 2>&1 | tee -a "$LOG" | grep -v "already satisfied"
  else
    "$VENV/bin/pip" install -r "$OPT_REQ" >>"$LOG" 2>&1
  fi
  if "$VENV/bin/python" -c "import celltypist, decoupler, liana" >>"$LOG" 2>&1; then
    echo "$OPT_SHA" >"$VENV/.optional.sha"
    log "Expert add-on packages installed."
  else
    # left unrecorded on purpose, so the next launch tries again; the analyses that need these
    # packages report themselves as skipped rather than failing
    log "WARN: expert add-on packages are not available — those analyses will report themselves as skipped."
  fi
fi

# ---- self-test: the server code must import cleanly with these packages
cd "$APP_DIR" || fail "App folder missing: $APP_DIR"
if ! "$VENV/bin/python" -c "import server.app" >>"$LOG" 2>&1; then
  fail "The app failed its start-up self-test (a package is incompatible). The log shows which one."
fi

# ---- optional: self-test instead of starting (Run self-test.command)
if [ "${1:-}" = "--selftest" ]; then
  echo ">> Running the self-test (about 2 minutes)…"
  "$VENV/bin/python" -m server.selftest 2>&1 | grep -v "Warning" | tee -a "$LOG"
  echo; read -r -p "Press Enter to close."
  exit 0
fi

# ---- start the server
log "Starting server on port $PORT"
if [ "$MODE" = terminal ]; then
  echo ">> Starting http://localhost:$PORT  (leave this window open; press Ctrl+C to stop)"
  ( for _ in $(seq 1 60); do curl -s -m 1 "http://localhost:$PORT/api/settings" >/dev/null 2>&1 && { open "http://localhost:$PORT"; exit 0; }; sleep 1; done ) &
  exec "$VENV/bin/python" -m uvicorn server.app:app --port "$PORT" 2>&1 | tee -a "$LOG"
fi

# GUI mode: the server runs on its own in the background; this launcher exits once it is up.
# Stop it from the sidebar ("Quit RNAseq Bench"), with "Stop RNAseq Bench.command", or just leave it — it uses
# no CPU while idle.
nohup "$VENV/bin/python" -m uvicorn server.app:app --port "$PORT" >>"$LOG" 2>&1 &
SERVER=$!
disown "$SERVER" 2>/dev/null
for _ in $(seq 1 90); do
  curl -s -m 1 "http://localhost:$PORT/api/settings" >/dev/null 2>&1 && break
  kill -0 "$SERVER" 2>/dev/null || break
  sleep 1
done
if ! curl -s -m 2 "http://localhost:$PORT/api/settings" >/dev/null 2>&1; then
  kill "$SERVER" 2>/dev/null
  fail "The server did not start. The last lines of the log say why."
fi
open "http://localhost:$PORT"
note "Running at http://localhost:$PORT"
exit 0
