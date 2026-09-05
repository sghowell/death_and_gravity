import sys
from pathlib import Path

PHYSICAL = Path(__file__).resolve().parents[1]
for source in (PHYSICAL/"src", PHYSICAL.parent/"src", PHYSICAL.parents[1]/"src"):
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
