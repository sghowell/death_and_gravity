import sys
from pathlib import Path

SCATTERING = Path(__file__).resolve().parents[1]
for source in (SCATTERING/"src", SCATTERING.parent/"control"/"src", SCATTERING.parent/"physical"/"src",
               SCATTERING.parent/"src", SCATTERING.parents[1]/"src"):
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
