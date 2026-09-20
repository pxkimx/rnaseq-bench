# Test harness

Synthetic datasets that mimic what GEO serves, and a runner that executes the pipelines exactly as the app does.

```
python tests/make_10x.py data10x            # 6 GSM-prefixed 10x v3 triplets + series matrix
python tests/make_more.py                   # RAW.tar, dense csv.gz, 10x .h5, processed .h5ad (from data10x/)
python tests/make_geo_variants.py           # nested per-sample tar.gz in RAW.tar, per-sample .h5, per-sample csv
python tests/make_bulk_variants.py          # per-sample HTSeq files in RAW.tar + series matrix, .xlsx + README
python tests/run_job.py auto data_nested job_nested          # run one job; prints flags / traceback
python tests/mock_geo_server.py mock_geo    # fake ftp.ncbi.nlm.nih.gov listing on :8799 (TL_GEO_BASE=http://127.0.0.1:8799)
python tests/fake_agent_test.py             # agent tool loop with a scripted fake Anthropic client
python -m server.selftest                   # what "Run self-test.command" runs
```
Tested with fresh virtualenvs on Python 3.10 and 3.12 (see requirements.txt header).
