import sys, os, json, types, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
os.environ.setdefault("TL_HOME", os.path.expanduser("~/Library/Application Support/RNAseqBench")); os.environ["ANTHROPIC_API_KEY"] = "sk-test"
import anthropic
from server import agent
# the scripted tool call reads analyzed.h5ad, so pick the most recent *single-cell* job
def _newest_sc():
    jobs = sorted((pathlib.Path(os.environ["TL_HOME"]) / "jobs").iterdir(), key=lambda p: p.stat().st_mtime, reverse=True)
    for d in jobs:
        if (d / "analyzed.h5ad").exists():
            return d.name
    raise SystemExit("no finished single-cell job in TL_HOME — run one first (tests/run_job.py sc ...)")


JOB = _newest_sc()
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
png = pathlib.Path(os.environ["TL_HOME"]) / "jobs" / JOB / "agent_umap.png"
print("agent_umap.png exists:", png.exists())
assert png.exists(), "run_python did not produce the figure"
png.unlink()
print(agent.available_models())
