"""Expose every frozen P8 source root without modifying the environment."""
import sys
from pathlib import Path

ROOT = next(parent for parent in Path(__file__).resolve().parents if (parent/"problems/P8").is_dir())
for source in sorted((ROOT/"problems/P8").rglob("src")):
    sys.path.insert(0, str(source))
