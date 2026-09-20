import sys, json, pathlib, shutil, time
sys.path.insert(0, "/home/claude/rnaseq-bench")
from server.common import Job, run_safely
from server import sc_pipeline, bulk_pipeline
from server.io_utils import guess_kind, unpack_archives
kind, src, dst = sys.argv[1], pathlib.Path(sys.argv[2]), pathlib.Path(sys.argv[3])
params = json.loads(sys.argv[4]) if len(sys.argv) > 4 else {}
if dst.exists(): shutil.rmtree(dst)
(dst / "input").mkdir(parents=True)
files = [shutil.copy(p, dst / "input" / p.name) for p in sorted(src.iterdir()) if p.is_file()]
files = unpack_archives([pathlib.Path(f) for f in files], dst / "input")
k = guess_kind(files) if kind == "auto" else kind
print("kind", k, "files", [f.name for f in files][:8])
job = Job(dst, k, "test")
(dst / "request.json").write_text(json.dumps({"kind": k, "name": "test", "params": params, "files": [str(f) for f in files]}))
t = time.time()
run_safely(job, sc_pipeline.run if k == "sc" else bulk_pipeline.run, files, params)
st = json.loads((dst / "status.json").read_text())
print(f"{time.time()-t:.0f}s", st["state"], st.get("error", "")); print(st.get("trace", "")[-3000:])
if st["state"] == "done":
    r = json.loads((dst / "result.json").read_text())
    print("flags:", *[f"[{f['level']}] {f['text'][:140]}" for f in r["flags"]], sep="\n  ")
    print("sections:", [(s["id"], len(s["items"])) for s in r["sections"]])
