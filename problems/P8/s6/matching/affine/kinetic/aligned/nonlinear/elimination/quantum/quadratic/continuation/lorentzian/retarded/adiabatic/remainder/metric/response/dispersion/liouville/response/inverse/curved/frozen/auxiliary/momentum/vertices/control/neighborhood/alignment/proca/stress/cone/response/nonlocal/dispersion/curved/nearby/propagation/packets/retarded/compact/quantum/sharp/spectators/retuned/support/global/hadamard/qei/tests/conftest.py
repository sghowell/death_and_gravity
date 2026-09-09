"""Make every frozen P8 source root available without installation."""
import sys
from pathlib import Path

ROOT=next(parent for parent in Path(__file__).resolve().parents if (parent/"problems/P8").is_dir())
for source in sorted((ROOT/"problems/P8").rglob("src")):
    sys.path.insert(0,str(source))
