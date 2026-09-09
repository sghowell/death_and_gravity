"""Native standalone import paths for this scalar and its immutable ancestry."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P8 = next(parent for parent in ROOT.parents if parent.name == "P8")
sys.path[:0] = [str(source) for source in sorted(P8.rglob("src"))]
