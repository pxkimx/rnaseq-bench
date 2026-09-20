#!/bin/bash
# Double-click this file: it opens Terminal so you can watch every step (useful if the app itself
# will not open). Same environment and data as the app.
cd "$(dirname "$0")" || exit 1
APP="$(pwd)/RNAseq Bench.app/Contents/Resources/app"
[ -d "$APP" ] || APP="$(pwd)"           # also works when placed next to server/ and web/
MODE=terminal exec /bin/bash "$APP/launch.sh"
