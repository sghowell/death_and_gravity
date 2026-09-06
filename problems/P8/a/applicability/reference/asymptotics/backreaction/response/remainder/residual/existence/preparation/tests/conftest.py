import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P8 = next(parent for parent in ROOT.parents if parent.name == "P8")
for source in P8.rglob("src"):
    if source.is_dir():
        sys.path.insert(0, str(source))
