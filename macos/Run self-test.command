#!/bin/bash
# Double-click: installs/updates the environment if needed, then runs a 2-minute simulated analysis
# (single-cell + bulk) to confirm everything works on this Mac. Nothing is started.
cd "$(dirname "$0")" || exit 1
APP="$(pwd)/RNAseq Bench.app/Contents/Resources/app"
[ -d "$APP" ] || APP="$(pwd)"
MODE=terminal exec /bin/bash "$APP/launch.sh" --selftest
