"""Read-only replay of S6.1 vacuum underdetermination, not UV admissibility."""

import argparse
import hashlib
import json
from pathlib import Path

from p8_scattering import uv_audit
from p8_scattering.verify import prior_checks as scattering_prior_checks

from . import dispersion, extension, vacuum

ROOT = Path(__file__).resolve().parents[2]
SCATTERING = ROOT.parent/"s5"/"scattering"
REPORT = ROOT/"certificates"/"vacuum-extension.json"
PRIOR_SHA = "086e5047f1751d4f28ddb6c813651b4b6d0aaf3471cf2ef9a7828ffcf038772f"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_prior():
    path = SCATTERING/"certificates"/"tree-control.json"
    if sha(path) != PRIOR_SHA:
        raise ValueError("Published S5.3b/4 checkpoint changed")
    prior = json.loads(path.read_text())
    for relative, expected in prior["source_sha256"].items():
        if sha(SCATTERING/relative) != expected:
            raise ValueError(f"Published S5.3b/4 input changed: {relative}")
    if scattering_prior_checks() != prior["prior_oscillator_sha256"]:
        raise ValueError("Published oscillator checkpoint changed")
    return prior


def build_report():
    prior = check_prior()
    audit = uv_audit.derive()
    if audit != prior["uv_applicability"]:
        raise ValueError("Vacuum obstruction differs from the pinned audit")
    sources = sorted(ROOT.glob("src/p8_uv/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S6.1.D", "date": "2026-09-04",
            "status": "VACUUM_UNDERDETERMINATION_VERIFIED; CONDITIONAL_UV_POLICY_ADOPTED; MATCHING_OPEN",
            "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
            "prior_scattering_sha256": PRIOR_SHA,
            "extension": extension.derive(), "vacuum_amplitude": vacuum.derive(),
            "original_polynomial_vacuum_audit": audit,
            "bound_bookkeeping_synthetic_controls_not_P8_bounds": {
                "missing_gravity": dispersion.necessary_bound(-1, error=0, gravity_lower_radius=None),
                "gravity_allows_negative": dispersion.necessary_bound(-1, error=0, gravity_lower_radius=2),
                "strict_negative_margin": dispersion.necessary_bound(-3, error=0, gravity_lower_radius=2)},
            "written_lemmas": ["notes/vacuum-extension.md"],
            "conditional_policy": "FORMULATION.md; assumptions are not machine-proved physics",
            "selected_followup": "S6.2: test a common matching construction for the positive-lambda member",
            "not_established": ["controlled matching between vacuum and clock tube",
                                "healthy transition-domain spectrum or heavy gap",
                                "finite-gravity or bounce positivity verdict", "loop/all-orders or radiative control",
                                "M1 interactions", "nonlinear/global stability", "UV completion", "P8(a)"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("S6.1 report differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.1: smooth vacuum freedom and conditional analytic obstruction replay passed; UV matching OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
