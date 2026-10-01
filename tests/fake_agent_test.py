import atexit, os, pathlib, shutil, sys, tempfile
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
os.environ["ANTHROPIC_API_KEY"] = "sk-test"

# The scripted tool call reads analyzed.h5ad, so this needs a finished single-cell analysis. It runs
# against a throw-away copy of one rather than against the real job: the agent writes agent_umap.png
# into the job folder, and even though the test removes it again, the folder's timestamp changes and
# the user's analysis jumps to the top of Recent with today's date on it. The files are hard-linked,
# so nothing is copied and the original's timestamps are not touched.
SRC_HOME = pathlib.Path(os.environ.get("RB_SRC_HOME", os.path.expanduser("~/Library/Application Support/RNAseqBench")))


def _newest_sc(home):
    d = home / "jobs"
    if not d.is_dir():
        raise SystemExit(f"no analyses in {home} — run one first (tests/run_job.py sc ...)")
    for j in sorted(d.iterdir(), key=lambda p: p.stat().st_mtime, reverse=True):
        if (j / "analyzed.h5ad").exists():
            return j
    raise SystemExit(f"no finished single-cell job in {home} — run one first (tests/run_job.py sc ...)")


src = _newest_sc(SRC_HOME)
TMP = pathlib.Path(tempfile.mkdtemp(prefix="rnaseq-bench-agent-test-"))
atexit.register(shutil.rmtree, TMP, True)
JOB = src.name
work = TMP / "jobs" / JOB
work.mkdir(parents=True)
for f in src.iterdir():
    if f.is_file():
        try:
            os.link(f, work / f.name)          # same filesystem: no copy, and the original is untouched
        except OSError:
            shutil.copy2(f, work / f.name)     # different filesystem: fall back to a copy
print(f"testing against a throw-away copy of {JOB} in {TMP}")
os.environ["TL_HOME"] = str(TMP)               # agent.JOBS is read at import, so this must come first

import anthropic                               # noqa: E402
from server import agent                       # noqa: E402
class Blk:
    def __init__(self, **k): self.__dict__.update(k)
    def model_dump(self): return dict(self.__dict__)
class Delta: 
    def __init__(self, t): self.type="text_delta"; self.text=t
class Ev:
    def __init__(self, t): self.type="content_block_delta"; self.delta=Delta(t)
class Stream:
    def __init__(self, msgs, n): self.n=n; self.msgs=msgs
    def __enter__(self): return self
    def __exit__(self,*a): pass
    def __iter__(self):
        if self.n == 0: yield Ev("Let me look. ")
        if self.n == 2: yield Ev("Done: 10 ISGs up.")
    def get_final_message(self):
        if self.n == 0: return Blk(content=[Blk(type="text", text="Let me look. "), Blk(type="tool_use", id="t1", name="get_result", input={"job": JOB})])
        if self.n == 1: return Blk(content=[Blk(type="tool_use", id="t2", name="run_python", input={"job": JOB, "code": "import scanpy as sc, matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt\na=sc.read_h5ad('analyzed.h5ad'); print(a.shape); sc.pl.umap(a, color='leiden', show=False); plt.savefig('agent_umap.png')"})])
        return Blk(content=[Blk(type="text", text="Done: 10 ISGs up.")])
class Msgs:
    n = 0
    def stream(self, **kw):
        s = Stream(kw["messages"], self.n); self.n += 1; return s
class Models:
    def list(self, limit=100): return [Blk(id="claude-sonnet-4-5"), Blk(id="claude-opus-4-1")]
class Fake:
    def __init__(self, api_key=None): self.messages=Msgs(); self.models=Models()
anthropic.Anthropic = Fake
out = list(agent.chat("c1", "what changed?", JOB, lambda j, p: {"job": "x"}))
for o in out: print(o.strip()[:160])
assert any('"type": "done"' in o for o in out)
png = work / "agent_umap.png"
print("agent_umap.png exists:", png.exists())
assert png.exists(), "run_python did not produce the figure"
assert not (src / "agent_umap.png").exists(), "the test wrote into the user's own analysis folder"
print(agent.available_models())
