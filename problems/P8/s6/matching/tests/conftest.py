import sys
from pathlib import Path

MATCHING = Path(__file__).resolve().parents[1]
P8 = MATCHING.parents[1]
for source in (MATCHING/"src", P8/"src", P8/"s5"/"src", P8/"s5"/"physical"/"src",
               P8/"s5"/"control"/"src", P8/"s5"/"scattering"/"src", P8/"s6"/"src"):
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
