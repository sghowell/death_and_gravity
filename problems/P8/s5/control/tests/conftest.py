import sys
from pathlib import Path

CONTROL = Path(__file__).resolve().parents[1]
for source in (CONTROL/"src", CONTROL.parent/"physical"/"src", CONTROL.parent/"src", CONTROL.parents[1]/"src"):
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
