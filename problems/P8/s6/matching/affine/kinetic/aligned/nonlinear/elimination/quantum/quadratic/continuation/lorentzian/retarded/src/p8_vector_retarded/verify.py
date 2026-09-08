"""Read-only selected-state constrained retarded mass-response checkpoint."""
import argparse
import json
from functools import cache
from pathlib import Path

import sympy as sp
from p8_affine import verify as affine
from p8_vector_lorentzian_pole import verify as parent

from . import bubble, canonical, causal, clock, flat, modes

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"selected-state-retarded-mass-response.json"
PARENT_SHA = "13944f63c00a27e4ef8d7d4f7d9e596fbdc8906e3b3014e42281e3a75a35da10"
sha, serialize = affine.sha, affine.serialize


def source_files():
    return (sorted(ROOT.glob("src/p8_vector_retarded/*.py"))+sorted(ROOT.glob("tests/*.py"))
            +sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md")))


@cache
def prior_checks():
    if sha(parent.REPORT) != PARENT_SHA:
        raise ValueError("The frozen physical-signature quadratic pole certificate changed")
    parent.validate_report(json.loads(parent.REPORT.read_text()), parent.build_report())
    return {"S6_63_fully_rebuilt": PARENT_SHA, "original_clock_state_and_fixed_profiles_unchanged": True}


@cache
def residuals():
    out = {}
    for module in (canonical, modes, bubble, flat, causal, clock):
        prefix = module.__name__.rsplit(".", 1)[1]
        out.update({prefix+"_"+name: value for name, value in module.checks().items()})
    return out


@cache
def proof_checks():
    data = clock.profiles()
    out = {"actual_temporal_mass_jet_positive": data["alpha"].is_positive,
           "actual_spatial_mass_jet_positive": data["beta"].is_positive,
           "selected_initial_surface_unchanged": data["initial_time"] == -sp.Rational(1, 2),
           "TT_angular_weight_positive": all(value >= 0 for value in sp.Poly(1+modes.cosine**2, modes.cosine).all_coeffs()),
           "mixed_angular_weight_is_real_square_on_unit_circle": modes.angular(1-modes.cosine**2-modes.sine**2) == 0,
           "three_physical_polarizations_only": all(len(modes.polarization(0, j)) == 3 for j in (0, 1)),
           "temporal_mode_determined_by_longitudinal_momentum":
               all(modes.polarization(0, j)[2][0].has(modes.data(0, j)["pL"]) for j in (0, 1))}
    control = flat.controls()
    point = {flat.alpha: 1, flat.beta: 0, flat.m: 1, flat.w: 2}
    out["omitted_contact_static_fixture_nonzero"] = control["omitted_longitudinal_contact_per_momentum"].subs(point) == sp.Rational(3, 4)
    out["omitted_contact_UV_fixture_nonzero"] = control["omitted_contact_UV_residue_times_pi_squared"].subs(point) == sp.Rational(3, 64)
    for channel in bubble.CHANNELS:
        out[channel+"_no_binary_float_in_kernel"] = not bubble.closed(channel).atoms(sp.Float)
    return {name: bool(value) for name, value in out.items()}


def controls():
    calls = []
    bad_labels = (True, False, sp.true, sp.false, 0.0, sp.Float(0), sp.Integer(0), "0", -1, 2, sp.Rational(1, 2))
    for value in bad_labels:
        calls.extend((lambda value=value: modes.data(value, 0),
                      lambda value=value: modes.data(0, value),
                      lambda value=value: modes.direction(value),
                      lambda value=value: modes.polarization(value, 0),
                      lambda value=value: modes.polarization(0, value),
                      lambda value=value: modes.wightman(value),
                      lambda value=value: bubble.vertex(value)))
    for value in (True, False, sp.true, 0, 1, "transverse", "longitudinal", "", None):
        calls.append(lambda value=value: canonical.insertion(value))
    for value in (True, False, sp.true, 0, 1, "T", "L", "XX", "", None):
        calls.extend((lambda value=value: bubble.direct(value),
                      lambda value=value: bubble.closed(value)))
    for value in (True, False, sp.true, sp.false, 0.0, sp.Float(0), sp.Integer(0), "0", -1, 3):
        calls.extend((lambda value=value: flat.coefficient(value),
                      lambda value=value: flat.integrated_residue(value)))
    for value in (True, False, sp.true, 1.0, sp.Float(1), sp.Integer(1), "1", -5, 0, 3):
        calls.append(lambda value=value: flat.radial_residue(value))
    calls.extend(lambda value=value: serialize(value) for value in (1.0, sp.Float(1), sp.oo, sp.nan))
    rejected = 0
    for call in calls:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(calls):
        raise ValueError("An unsupported retarded mass-response input was accepted")
    return {"rejected_inputs": rejected, "native_index_validation_before_cache": True,
            "temporal_constraint_contact_retained": True,
            "mass_only_source_not_full_metric_lapse": True,
            "finite_cutoff_not_claimed_microcausal": True,
            "retarded_background_support_not_finite_perturbed_cone": True,
            "no_early_dimensional_finite_limit": True}


@cache
def build_report():
    prior, identities = prior_checks(), residuals()
    exact, proofs = affine.certify_residuals(identities), proof_checks()
    if not all(value is True for value in proofs.values()):
        raise ValueError("A constrained retarded mass-response audit check failed")
    scalar_count = sum(value.rows*value.cols if isinstance(value, sp.MatrixBase) else 1 for value in identities.values())
    return {"schema": 1, "claim": "P8-S6.64.RETARDED_MASS_RESPONSE", "date": "2026-09-08",
            "status": "BARE_SELECTED_STATE_MASS_RESPONSE_AND_CONTACT; ORIGINAL_P8_OPEN",
            "prior_sha256": prior, "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in source_files()},
            "formulation": "FORMULATION.md", "written_proofs": ["notes/constraint.md", "notes/bubble.md", "notes/causal.md"],
            "exact_residuals": exact, "named_exact_check_count": len(identities), "checked_scalar_entries": scalar_count,
            "proof_checks": proofs,
            "unchanged_clock_profiles": serialize(clock.profiles()),
            "constraint_and_contact": serialize(canonical.constraint()),
            "canonical_insertions": serialize({sector: canonical.insertion(sector) for sector in ("T", "L")}),
            "general_momentum_connected_integrands": serialize({sector: bubble.closed(sector) for sector in bubble.CHANNELS}),
            "mode_normalization": "Displayed U vectors in code are a^(3/2) times the physical orthonormal modes. V=m^2 diag(alpha,-beta,-beta,-beta). Pair vertices include a^3 and the density-volume factors cancel. Two transverse polarizations and one constrained longitudinal polarization are retained.",
            "retarded_formula": "For this parity-invariant selected state, R_bubble=-2 theta(u-s) Im integral d^3p/(2pi)^3 sum C_sector^+(u,s,p,K-p). Add -delta(u-s) integral d^3p/(2pi)^3 alpha(u)^2 z(u,p)|p_L(u,p)|^2, independent of external K. These are coordinate-density responses; normalize the output by a(u)^3 when comparing to the physical normalized counterterm.",
            "canonical_covariance_identity": "W=F(Sigma+iJ/2), C^+=Tr(A W B W^T)/2, delta M=-J B. i(C^+-conj(C^+))=Tr(A F(delta M Sigma+Sigma delta M^T)F^T)/2. All real symmetric input/output/covariance matrices and general real F are symbolic; no stationarity is assumed.",
            "flat_frequency_kernel": serialize(flat.frequency_kernel()),
            "flat_pole_and_omitted_contact_controls": serialize(flat.controls()),
            "prepared_source_and_causal_boundary": "n is smooth and compactly supported strictly after u0=-1/2, zero on an initial neighborhood. This fixes the original initial covariance for this restricted perturbative mass-source response only. Ordinary-Proca CCR prove background-cone support away from coincidence. Local contacts/counterterms do not enlarge that support; they do not establish a finite off-clock cone or an arbitrary compatible varied-state family.",
            "controls": controls(),
            "verdict": "The constrained selected-state bare mass-response integrands, temporal contact, exact homogeneous covariance identity and all three flat time-frequency UV poles agree. The finite matched causal kernel and original P8 remain open.",
            "not_established": ["The dimensionally regulated combined finite bare/contact/local response, finite matching or integrated feedback bounds",
                                "Full metric and second mass-source blocks, general varied Cauchy covariance, or finite off-clock cones",
                                "Coupled quantum stability, interacting cutoff, common UV parent, V/G/B or original P8 closure"],
            "verification_boundary": "Exact canonical, polarization, Wick, covariance, Gamma-residue and frozen Feynman comparisons, with written constrained-field and causal-support arguments; not proof-assistant formalization"}


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("The constrained retarded mass-response report differs from read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.64.RETARDED_MASS_RESPONSE replay passed; finite kernel and original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
