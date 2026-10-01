"""A job's status.json must never be caught half written.

status.json is rewritten at every step of an analysis while the page polls it every 1.5 s and the
assistant reads it directly. A plain write_text truncates the file before writing it, so a reader can
open it empty: in the sibling tool MassSpec Bench that surfaced as a 500 from the status endpoint,
mid-analysis, for no reason the user could see.

Run:  python3 tests/test_job_files.py
"""
import json
import shutil
import sys
import tempfile
import threading
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from server.common import write_json  # noqa: E402

fails = []


def check(name, got, want):
    ok = got == want
    if not ok:
        fails.append(f"{name}: got {got!r}, wanted {want!r}")
    print(f"{'ok  ' if ok else 'FAIL'} {name}")


def hammer(path: Path, writer, seconds=2.0):
    """One thread rewrites the file as fast as it can, another reads it. Returns (reads, torn reads)."""
    stop = threading.Event()
    counts = {"reads": 0, "torn": 0}

    def write():
        i = 0
        while not stop.is_set():
            i += 1
            writer(path, {"state": "running", "pct": i % 100, "message": "Step " + "x" * (i % 400),
                          "time": time.time()})

    def read():
        while not stop.is_set():
            try:
                json.loads(path.read_text())
            except FileNotFoundError:
                counts["torn"] += 1          # the file vanished mid-swap — also a torn read
            except json.JSONDecodeError:
                counts["torn"] += 1
            counts["reads"] += 1

    writer(path, {"state": "running", "pct": 0, "message": "start", "time": time.time()})
    w = threading.Thread(target=write, daemon=True)
    r = threading.Thread(target=read, daemon=True)
    w.start(); r.start()
    time.sleep(seconds)
    stop.set(); w.join(1); r.join(1)
    return counts["reads"], counts["torn"]


tmp = Path(tempfile.mkdtemp(prefix="rnaseq-bench-jobfiles-"))
try:
    # the control: what the code used to do, so this test fails loudly if the helper is ever
    # quietly replaced by a plain write again
    def plain(path, data):
        path.write_text(json.dumps(data))

    reads, torn = hammer(tmp / "plain.json", plain)
    print(f"     plain write_text:  {torn} torn out of {reads} reads")
    check("a plain write really can be caught half written", torn > 0, True)

    reads, torn = hammer(tmp / "status.json", write_json)
    print(f"     write_json:        {torn} torn out of {reads} reads")
    check("write_json is never caught half written", torn, 0)
    check("the reader actually got to read", reads > 100, True)

    # and the file left behind is valid, with no .tmp beside it
    final = json.loads((tmp / "status.json").read_text())
    check("the file left behind is complete", final["state"], "running")
    check("no temporary file is left behind", sorted(p.name for p in tmp.iterdir()), ["plain.json", "status.json"])
finally:
    shutil.rmtree(tmp, ignore_errors=True)

print()
if fails:
    print(f"{len(fails)} failed:")
    for f in fails:
        print("  " + f)
    sys.exit(1)
print("all job-file checks passed")
