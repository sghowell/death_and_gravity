"""Read-only reciprocal-scale source, geometry and affine-bound replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

from p8_variable_beta import verify as prior

from . import identities

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"reciprocal-affine-obstruction.json"
PRIOR_SHA = "335cd52028baf56b30c75377aa7e6edd7ec86db2799f49b13f467242674fc4d0"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


@cache
def prior_checks():
    if sha(prior.REPORT) != PRIOR_SHA:
        raise ValueError("The frozen source-aware variable-beta ancestor changed")
    prior.validate_report(json.loads(prior.REPORT.read_text()), prior.build_report())
    return {"S6_20_source_aware_metric_and_shift_variations_with_recursive_pins": PRIOR_SHA}


def controls():
    rejected = 0
    for index in range(4):
        for bad in (0, -1, True, 0.01):
            values = [2, 1, 1, 1]
            values[index] = bad
            try:
                identities.affine_tail_bound(*values)
            except (TypeError, ValueError):
                rejected += 1
    if rejected != 16:
        raise ValueError("An invalid tail-bound calibration input bypassed validation")
    return {"invalid_exact_bound_inputs_rejected": rejected,
            "unit_momentum_bound_example": str(identities.affine_tail_bound(2, Fraction(5, 4), Fraction(3, 2))),
            "momentum_three_bound_example": str(identities.affine_tail_bound(2, Fraction(5, 4), Fraction(3, 2), 3)),
            "numeric_boundary": "Positive inputs calibrate a conditional bound; they do not certify a background solution"}


@cache
def build_report():
    previous = prior_checks()
    residuals = identities.checks()
    audit = ROOT/"tests"/"test_reciprocal_independent_audit.py"
    if not audit.is_file():
        raise ValueError("The separately authored source/geometric audit must exist before freezing")
    sources = sorted(ROOT.glob("src/p8_reciprocal_geometry/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S6.25.RECIPROCAL", "date": "2026-09-07",
        "status": "CONDITIONAL_RECIPROCAL_F_AFFINE_INCOMPLETENESS_WITH_PHYSICAL_G_COMPLETENESS; ORIGINAL_P8_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_source_and_differential_residuals": dict.fromkeys(residuals, "0"),
        "exact_residual_count": len(residuals),
        "hypotheses": {
            "geometry": "Positive C2 a; h=a'/a; h0=0,h'0>0; (a*a')'>=0; b=alpha/a,alpha>0",
            "lapse_and_link": "Positive C1 c, finite at every finite u; positive continuous finite P=2(beta1+2beta2*y+beta3*y²)",
            "source": "No additional f null stress; actual undivided f metric equation retained",
            "normalization": "Dimensionless g time u=T/tau; beta=tau² beta_physical/M²; common positive Einstein coefficient divided out",
            "domain": "Specified future coordinate domain u>=0 with spatial R³; time-reflected hypotheses for past claims",
        },
        "continuous_proof_chain": {
            "center": "z=h/c, z0=0, z'0=h'0/c0>0; actual source equation forces c0>y0 without dividing byh",
            "invariant": "D=h/y-z; D'+[P/(2y²z)]D=(h/y)'>=0, integrated only from a positive puncture",
            "order": "D>0, c>y, z'>0 for every u>0; no finite future branch crossing",
            "scale": "(a*a')'>=0 and nondegenerate bounce imply a>=a0 and a tends to infinity",
            "finite_window": "Integral_u0^L b*c du <=alpha/z(u0)*(1/a(u0)-1/a(L))",
            "tail": "Integral_u0^infinity b*c du <=alpha/(z(u0)*a(u0)), for every u0>0",
            "geodesics": "f null dλ/du=bc/|p| and noncomoving timelike ds/du<=bc/|p|; g causal integrands bounded below positively",
            "CD_geometry": "For a=(1+u²)²,alpha2: (h/y)'=2(1+u²)²(1+7u²)>0",
        },
        "checked_controls": controls(),
        "not_established": [
            "No requirement that original P8 make the second metric complete; physical g is complete here",
            "No curvature singularity or common extension/inextendibility verdict",
            "No conclusion about comoving f timelike completeness or global scalar health",
            "No assertion for nonreciprocal scales, extra f null stress, P zeros/sign changes or dropped geometric hypothesis",
            "No original C/D operator map, light-only exclusion, mass gap, cutoff, loop or UV verdict",
            "No G1 certificate dependency or transfer of its separate mass asymptotics",
        ],
        "verification_boundary": "Exact source and differential identities plus separately authored audit and continuous written integrating-factor/geodesic proof; not proof-assistant formalized",
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Reciprocal geometric certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.25.RECIPROCAL: source/geometric replay passed; physical-g bounce and original P8 remain unexcluded")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
