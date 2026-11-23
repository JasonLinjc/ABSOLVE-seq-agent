"""Local inputs and optional caches. All overrides accept absolute paths."""
import os
import shutil
from pathlib import Path


def use_local(executor):
    if executor == "local":
        return True
    return not os.environ.get("SLURM_JOB_ID") and shutil.which("sbatch") is None


def location(variable, default):
    return Path(os.environ.get(variable, default)).expanduser().resolve()


ROOT = Path(__file__).resolve().parents[1]
RUNS = location('ABSOLVE_RUNS', ROOT / 'runs')
DATA = location('ABSOLVE_DATA', ROOT / 'data')
REFERENCES = DATA / 'references'
PUBLIC = REFERENCES / 'target_info'
MANIFEST = location('ABSOLVE_MANIFEST', DATA / 'sample_manifest.tsv')
DEPLOYMENT = location('ABSOLVE_CACHE', DATA / 'cache')
PUBLISHED = REFERENCES / 'published'
RUNTIME = ROOT / 'runtime'
