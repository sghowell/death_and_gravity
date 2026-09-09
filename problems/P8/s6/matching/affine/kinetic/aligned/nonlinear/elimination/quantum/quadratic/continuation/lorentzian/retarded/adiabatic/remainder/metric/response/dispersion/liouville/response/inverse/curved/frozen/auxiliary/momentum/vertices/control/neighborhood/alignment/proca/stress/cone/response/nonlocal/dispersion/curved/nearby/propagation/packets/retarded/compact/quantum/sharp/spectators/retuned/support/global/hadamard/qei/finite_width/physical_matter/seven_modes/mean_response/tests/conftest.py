"""Make immutable P8 source roots available without installation."""

import sys
from pathlib import Path

ROOT = next(p for p in Path(__file__).resolve().parents if (p / "problems/P8").is_dir())
for source in sorted((ROOT / "problems/P8").rglob("src")):
    sys.path.insert(0, str(source))
