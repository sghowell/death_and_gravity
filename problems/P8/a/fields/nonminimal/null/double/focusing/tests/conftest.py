"""Native standalone imports without changing ancestor scientific code."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
P8 = next(p for p in ROOT.parents if p.name == "P8")
sys.path[:0] = [str(p) for p in sorted(P8.rglob("src"))]
