"""Separate import roots; do not change the pinned earlier checkpoints."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT.parent / "src"))
sys.path.insert(0, str(ROOT.parents[1] / "src"))
