#!/bin/bash
# Build / refresh the Mac app from this source tree.
#   ./build_mac.sh            -> updates the app bundle in place (default: ../RNAseq Bench.app when it exists,
#                                otherwise macos/RNAseq Bench.app inside this folder) and writes a versioned zip
#   OUT_DIR=~/Desktop ./build_mac.sh   -> put the zip somewhere else
# After building, double-click the app (or run launch.sh) — the launcher notices the new version and swaps it in.
set -e
cd "$(dirname "$0")"
V=$(cat VERSION)
if [ -n "${APP:-}" ]; then :; elif [ -d "../RNAseq Bench.app" ]; then APP="../RNAseq Bench.app"; else APP="macos/RNAseq Bench.app"; fi
OUT_DIR="${OUT_DIR:-$(dirname "$APP")}"
[ -d "$APP/Contents" ] || { echo "no app bundle at $APP — copy macos/RNAseq Bench.app there first"; exit 1; }
sed -i.bak "s|<key>CFBundleVersion</key><string>[^<]*</string>|<key>CFBundleVersion</key><string>$V</string>|; s|<key>CFBundleShortVersionString</key><string>[^<]*</string>|<key>CFBundleShortVersionString</key><string>$V</string>|" "$APP/Contents/Info.plist" && rm -f "$APP/Contents/Info.plist.bak"
mkdir -p "$APP/Contents/Resources/app"
rm -rf "$APP/Contents/Resources/app/server" "$APP/Contents/Resources/app/web"
cp -R server web requirements.txt README.md VERSION CHANGELOG.md launch.sh "$APP/Contents/Resources/app/"
find "$APP" -name __pycache__ -type d -exec rm -rf {} + 2>/dev/null || true
chmod +x "$APP/Contents/MacOS/"* "$APP/Contents/Resources/app/launch.sh"
echo "app updated: $APP (v$V)"
if command -v zip >/dev/null 2>&1; then
  OUT="$OUT_DIR/RNAseq Bench v$V.zip"; rm -f "$OUT"
  ( cd "$(dirname "$APP")"
    items=("$(basename "$APP")")
    for f in *.command README-FIRST.txt; do [ -e "$f" ] && items+=("$f"); done
    zip -qr "$OUT" "${items[@]}" )
  echo "zip: $OUT"
fi
