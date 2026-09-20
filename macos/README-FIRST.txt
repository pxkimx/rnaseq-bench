RNAseq Bench — how to start
==============================

1. Unzip, then move "RNAseq Bench.app" anywhere (Applications or Desktop).

2. First open only — macOS Gatekeeper will say the app is "damaged" or "cannot be verified"
   because it is not signed with a paid Apple developer certificate. Fix it once:
      open Terminal, type   xattr -cr    (with a space after it), drag the app into the
      Terminal window so its path is filled in, press Return.
   Then double-click the app. (Alternative: right-click the app -> Open -> Open.)

3. The first launch installs the analysis packages into
   ~/Library/Application Support/RNAseqBench   (3–5 minutes, needs internet).
   Later launches take a few seconds. Your browser opens http://localhost:8765 and the
   server keeps running quietly in the background (no window). To stop it: "Quit RNAseq
   Bench" in the app's sidebar, or double-click "Stop RNAseq Bench.command".

Updating to a new version: just replace the .app (and run xattr -cr on it again). The
launcher notices the version change, stops any older copy that is still running, and
installs any new packages automatically. Your analyses and API key are kept.

If the app will not open at all: double-click "Start RNAseq Bench.command" instead —
it runs the same thing inside a Terminal window so you can see every message.
The log is at ~/Library/Application Support/RNAseqBench/rnaseq-bench.log

Needs Python 3.10–3.13 (python.org installer recommended; 3.14 is not supported by the analysis libraries yet — it can be installed alongside). The app tells you if it is missing.
