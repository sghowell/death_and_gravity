"""Read-only S5.3a exact oscillator/adiabatic-indicator replay, not a cutoff."""

import argparse
import hashlib
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_physical.verify import prior_checks

from . import oscillator as o

ROOT = Path(__file__).resolve().parents[2]
PHYSICAL = ROOT.parent/"physical"
REPORT = ROOT/"certificates"/"adiabatic.json"
PRIOR_SHA = "fec2b468d9144817026e923eeed1d7b573fdae068213b42a93cb2e4e77eb8c20"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def check_inputs():
    report = PHYSICAL/"certificates"/"interactions.json"
    if sha(report) != PRIOR_SHA:
        raise ValueError("Physical interaction checkpoint changed")
    for relative, expected in json.loads(report.read_text())["source_sha256"].items():
        if sha(PHYSICAL/relative) != expected:
            raise ValueError(f"Physical interaction input changed: {relative}")
    prior_checks()
    return PRIOR_SHA


@cache
def build_report():
    prior = check_inputs()
    residuals = o.generic_checks()
    if any(sp.cancel(value) != 0 for value in residuals.values()):
        raise ValueError("Nonzero oscillator normal-form residual")
    data = {}
    for chart in ("unitary", "gamma", "tensor"):
        out = o.derive(chart)
        data[chart] = {key: str(value) for key, value in out.items()}
    sources = sorted(ROOT.glob("src/p8_control/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {"schema": 1, "claim": "P8-S5.3a.D", "date": "2026-09-04",
            "status": "LINEAR_ADIABATIC_INDICATORS_BOUNDED; INTERACTION_CONTROL_OPEN",
            "prior_physical_sha256": prior,
            "source_sha256": {str(p.relative_to(ROOT)): sha(p) for p in sources},
            "exact_residuals": {key: "0" for key in residuals},
            "exact_volume_normalized_oscillators": data,
            "all_time_bounds": o.bounds(),
            "normal_form": "Y=a^(3/2)*sqrt(2*K)*Q, L=Ydot^2/2-Omega^2*Y^2/2 up to a time boundary",
            "local_units": "Omega^2*ell^2=q+mu, q=ell^2*k_physical^2; ell constant per local patch",
            "domain": "unitary x^2>=1/65, gamma x^2<=1/17 with q>=68, tensor all x; q is local, not fixed comoving",
            "written_lemmas": ["notes/oscillators.md"],
            "not_established": ["a global-in-time WKB solution error or vacuum prescription",
                                "interaction-coefficient adiabatic errors", "cubic-exchange or 2-to-2 amplitudes",
                                "a weak-coupling or cutoff hierarchy", "zero-mode/forward-channel treatment",
                                "all-orders/loops", "M1", "nonlinear stability", "UV completion"]}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("S5.3a certificate differs from exact replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S5.3a: exact oscillators and uniform adiabatic-indicator bounds replay passed; interaction control OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
