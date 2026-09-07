"""Standalone import paths for the unchanged A18/A16 certificate lineage."""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MAXWELL = ROOT.parent
FOCUSING = MAXWELL/"focusing"
COSMOLOGY = FOCUSING/"cosmology"
CLOCK = MAXWELL.parents[1]/"applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/qsei"
sys.path[:0] = [str(path/"src") for path in (ROOT, MAXWELL, FOCUSING, COSMOLOGY,
                                            CLOCK.parents[2]/"qsei", CLOCK, *CLOCK.parents)]
