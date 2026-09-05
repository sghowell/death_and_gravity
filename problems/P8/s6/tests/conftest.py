import sys
from pathlib import Path

S6 = Path(__file__).resolve().parents[1]
P8 = S6.parent
for source in (S6/"src", P8/"src", P8/"s5"/"src", P8/"s5"/"physical"/"src",
               P8/"s5"/"control"/"src", P8/"s5"/"scattering"/"src"):
    if str(source) not in sys.path:
        sys.path.insert(0, str(source))
