"""Read-only finite HR-tree sign theorem, exact controls and ancestor replay."""

import argparse
import hashlib
import json
from fractions import Fraction
from functools import cache
from pathlib import Path

import sympy as sp
from p8_bimetric_monotonic import verify as bimetric_prior
from p8_star import verify as star_prior

from . import components, controls, edge, graph, independent

ROOT = Path(__file__).resolve().parents[2]
REPORT = ROOT/"certificates"/"hr-tree-no-bounce.json"
P8 = next(path for path in ROOT.parents if path.name == "P8")
STAR_SHA = "54f23dfd05d7580969550d06ac4a72adc5faa61988be9e7c95eeb199fb546353"
BIMETRIC_SHA = "0efd16fa4f02f45f35448e2056fb37c863c29c2ed09203c3351ba9515baff031"
CONTRACT_SHA = "d00df35d5821b05baa61e897725a22ed1ac21d0536b74bdaece3ab0093215901"
CONVENTION_SHA = "73dd18f7a4b56a0205f9fbc6a4b09b213a3c280b3f0274414a255635b04e2b09"
CD_SHA = "caf8c8e688a7565b9d00f921c099a28da00f97522ed26ad182a8227eb80cd4dd"


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def serialize(value):
    """Canonical JSON data: no rounded arithmetic or floating enclosures."""
    if isinstance(value, (Fraction, sp.Basic)):
        return str(value)
    if isinstance(value, dict):
        return {str(key): serialize(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [serialize(item) for item in value]
    if value is None or isinstance(value, (str, int, bool)):
        return value
    raise TypeError("Only exact report values are allowed")


@cache
def prior_checks():
    pins = {"S6_16_constant_star_action_and_controls": (star_prior.REPORT, STAR_SHA),
            "S6_5_bimetric_monotonicity_context": (bimetric_prior.REPORT, BIMETRIC_SHA),
            "adopted_S6_matching_contract": (P8/"s6"/"FORMULATION.md", CONTRACT_SHA),
            "original_P8_B_curvature_dictionary": (P8/"FORMULATION.md", CONVENTION_SHA),
            "original_CD_matter_witness_context": (P8/"certificates"/"witness-CD_matter.json", CD_SHA)}
    if any(sha(path) != expected for path, expected in pins.values()):
        raise ValueError("A pinned action, convention, target or matching contract changed")
    star_prior.validate_report(json.loads(star_prior.REPORT.read_text()), star_prior.build_report())
    bimetric_prior.validate_report(json.loads(bimetric_prior.REPORT.read_text()), bimetric_prior.build_report())
    return {key: expected for key, (_, expected) in pins.items()}


def independent_bridges():
    stresses = []
    stress_keys = ("P", "rho_i", "pressure_i", "rho_j", "pressure_j", "null_i", "null_j")
    for fixture in independent.stress_fixtures():
        actual = edge.evaluate(fixture["beta"], fixture["y"], fixture["c"])
        if any(actual[key] != fixture[key] for key in stress_keys):
            raise ValueError("Independent four-coframe derivative versus primary stress failed")
        stresses.append({"y": str(fixture["y"]), "c": str(fixture["c"]),
                         "primary_minus_independent": dict.fromkeys(stress_keys, "0")})
    incidence = []
    for fixture in independent.incidence_fixtures():
        actual = graph.solve_fluxes(fixture["vertices"], fixture["edges"], fixture["divergences"])
        if actual != fixture["fluxes"]:
            raise ValueError("Independent incidence Gaussian elimination versus leaf solve failed")
        incidence.append({"vertices": fixture["vertices"],
                          "primary_minus_independent_fluxes": ["0"]*len(actual)})
    clocks = []
    clock_keys = ("A", "A_fixed_prime", "weighted_source", "lambda", "eta", "Hprime", "component_log_rate_bound", "yprimes")
    for fixture in independent.component_fixtures():
        actual = components.reconstruct(fixture["G"], fixture["y"], fixture["c"], fixture["nulls"],
                                        fixture["H"], fixture["component"], fixture["root"])
        if any(actual[key] != fixture[key] for key in clock_keys):
            raise ValueError("Independent component first-jet clock/source reconstruction failed")
        clocks.append({"vertices": len(fixture["G"]), "root": fixture["root"], "component": list(fixture["component"]),
                       "primary_minus_independent": dict.fromkeys(clock_keys, "0")})
    return {"literal_edge_stresses": stresses, "Gaussian_vs_leaf_incidence": incidence,
            "fixed_component_first_jets": clocks}


def checked_controls():
    d = edge.derive()
    point = dict(zip(d["beta"], (0, 1, 0, 0, 0), strict=True))
    point.update({d["y"]: 2, d["c"]: 3, d["H_i"]: 1, d["H_j"]: 0})
    wrong = {key: sp.factor(value.subs(point)) for key, value in edge.controls().items()}
    if tuple(wrong.values()) != (-12, 4, 1):
        raise ValueError("A reciprocal lapse-weight or factor-two omission failed to fire")
    mixed = controls.mixed_algebraic_tree()
    if sp.simplify(mixed["false_all_vertex_null_residual"]-2/(3*mixed["T"]**2)) != 0:
        raise ValueError("The actual algebraic-leaf all-vertex-inertia control failed")
    active = graph.active_component(3, mixed["edges"], mixed["P_values"])
    if active != (0, 2) or mixed["correct_component_null_residual"] != 0:
        raise ValueError("The actual selected-component control failed")
    if controls.zero_target_exception()["H_root_prime_at_zero"] != 2:
        raise ValueError("The actual zero-target Einstein exception failed")
    if controls.nec_violating_control()["null_at_bounce"] != -4:
        raise ValueError("The conserved non-NEC bounce control failed")
    cycle = ((0, 1), (1, 2), (2, 0))
    cycle_divergence = [sum((int(v == i)-int(v == j)) for i, j in cycle) for v in range(3)]
    if cycle_divergence != [0, 0, 0]:
        raise ValueError("The unit circulating cycle-flux countercontrol failed")
    nonconserved = graph.solve_fluxes(3, [(0, 1), (2, 1)], [1, -1, 0])
    if nonconserved != (1, 0):
        raise ValueError("Nonzero vertex-divergence omission control failed")
    tuned = edge.evaluate((-3, 1, 0, 0, -1), 1, 1)
    if tuned["P"] != 2 or any(tuned[key] for key in ("rho_i", "pressure_i", "rho_j", "pressure_j")):
        raise ValueError("The exact connected static endpoint-equality witness failed")
    endpoint = components.cd_endpoint_test(sp.Rational(1, 2), 0)
    equality = components.cd_endpoint_test(sp.Rational(1, 2), sp.Rational(8, 5))
    if endpoint["necessary_normalized_endpoint_error"] != sp.Rational(8, 5) or not endpoint["status"].startswith("EXCLUDED") or equality["status"] != "INCONCLUSIVE":
        raise ValueError("The strict actual-clock CD endpoint margin failed")
    invalid = [lambda: graph.validate_tree(3, cycle),
               lambda: graph.validate_tree(3, [(0, 1)]),
               lambda: graph.validate_tree(2, [(0, 1), (1, 0)]),
               lambda: graph.validate_tree(True, []),
               lambda: graph.active_component(2, [(0, 1)], [1], root=True),
               lambda: graph.solve_fluxes(2, [(0, 1)], [1, 0]),
               lambda: edge.evaluate((0, 1.0, 0, 0, 0), 1, 1),
               lambda: edge.evaluate((0, 1, 0, 0, 0), 1, 0),
               lambda: components.reconstruct([0, 1], [1, 1], [1, 1], [0, 0], 0, [0]),
               lambda: components.reconstruct([1, -1], [1, 1], [1, 1], [0, 0], 0, [0, 1]),
               lambda: components.reconstruct([1], [1], [1], [-1], 0, [0]),
               lambda: components.reconstruct([1], [1], [1], [0], True, [0]),
               lambda: components.reconstruct([1], [1], [1], [0], sp.oo, [0]),
               lambda: components.reconstruct([1], [1], [1], [0], 0, [0, 0]),
               lambda: components.cd_endpoint_test(sp.Rational(1, 2), -1)]
    rejected = 0
    for call in invalid:
        try:
            call()
        except (TypeError, ValueError):
            rejected += 1
    if rejected != len(invalid):
        raise ValueError("An invalid theorem-domain input was admitted")
    return serialize({"wrong_weight_and_normalization_defects": wrong,
                      "actual_algebraic_active_component": active,
                      "actual_false_all_vertex_inertia_residual": mixed["false_all_vertex_null_residual"],
                      "actual_zero_target_Hprime_at_zero": 2, "actual_nonNEC_null_at_bounce": -4,
                      "cycle_unit_flux_divergences_algebraic_only": cycle_divergence,
                      "nonconserved_vertex_fluxes_algebraic_only": nonconserved,
                      "static_tuned_connected_tree_edge_stresses": tuned,
                      "strict_CD_endpoint_test": endpoint, "CD_endpoint_equality_control": equality,
                      "invalid_domain_calls_rejected": rejected})


@cache
def build_report():
    previous = prior_checks()
    residuals = {"literal_edge_variations_and_reciprocal_weights": edge.checks(),
                 "literal_vertex_EH_sources_and_zero_EH_identity": edge.vertex_checks(),
                 "singleton_component": components.checks(1),
                 "three_vertex_fixed_component": components.checks(3),
                 "actual_background_and_separate_source_controls": controls.checks()}
    if any(sp.simplify(value) != 0 for group in residuals.values() for value in group.values()):
        raise ValueError("A finite-tree action, clock, component or actual-solution identity failed")
    if components.derive(3)["A"].is_positive is not True or components.derive(3)["eta"].is_nonnegative is not True:
        raise ValueError("Positive target inertia or nonnegative weighted NEC source failed")
    reserved = ROOT/"tests"/"test_hr_tree_independent_audit.py"
    if not reserved.is_file():
        raise ValueError("The separately authored covariant audit must be present before freezing")
    sources = sorted(ROOT.glob("src/p8_hr_tree/*.py"))+sorted(ROOT.glob("tests/*.py"))
    sources += sorted(ROOT.glob("*.md"))+sorted(ROOT.glob("notes/*.md"))
    return {
        "schema": 1, "claim": "P8-S6.17.TREE", "date": "2026-09-06",
        "status": "FINITE_HR_TREE_SEPARATE_NEC_PHYSICAL_SIGN_PRESERVATION; ORIGINAL_P8_OPEN",
        "prior_context_sha256": previous,
        "source_sha256": {str(path.relative_to(ROOT)): sha(path) for path in sources},
        "formulation": "FORMULATION.md", "written_proof": "notes/proof.md", "source_audit": "notes/sources.md",
        "exact_residuals": {name: dict.fromkeys(group, "0") for name, group in residuals.items()},
        "independent_Fraction_replay": serialize(independent.checks()),
        "primary_independent_bridges": independent_bridges(),
        "oriented_tree_pruning_fixtures": serialize(graph.checks()), "checked_controls": checked_controls(),
        "action_and_domain": {
            "interaction": "-2 sum_oriented_edges sqrt|g_i| sum_n beta_en e_n(sqrt(g_i^-1 g_j)); finite constant real beta0..4 of mass dimension four",
            "Einstein_coefficients": "Every constant G_v>=0; distinguished actual physical vertex has G_r>0; zero-EH vertices retain all their algebraic equations",
            "signature_and_curvature": "+---; R_B=-6(DH+2H²); EH=-G R_B/2; 3G H²=rho, -2G DH=rho+p",
            "old_normalization_dictionary": "S6.16 uses -2 beta; S6.5 m4*beta_old=2*beta_here",
            "geometry": "Finite connected acyclic tree; all actual metrics smooth regular common spatially flat FLRW with positive finite lapses/scales on a connected interval",
            "physical_metric_and_clock": "The actual vertex metric g_r; dT=N_r dt; component y_i=a_i/a_r,c_i=N_i/N_r, NOT c_i/y_i",
            "sources": "Each independent matter action couples to one vertex only; separately conserved on shell, isotropic and n_i=rho_i+p_i>=0; absent sources zero",
        },
        "literal_identities": {
            "edge_P": "P=2(beta1+2beta2*y+beta3*y²)",
            "edge_nulls": "I_i=(y-c)P; I_j=(c-y)P/(c*y³)",
            "one_lapse_null_reciprocity": "N_i*a_i³*I_i+N_j*a_j³*I_j=0",
            "source_aware_balance": "C_i=D_i rho_ie+3H_i I_i=3Pc(yH_j-H_i); C_j=-C_i/(c²y³)",
            "two_lapse_Bianchi_flux": "F_e=N_i²*a_i³*C_i=3N_i*N_j*a_i²*P*(a_jH_j-a_iH_i)=-N_j²*a_j³*C_j",
            "tree_incidence": "B_tree F=0 and finite-tree injectivity imply F_e=0 for every edge",
            "undivided_lock": "P_e(y_e)*(a_jH_j-a_iH_i)=0; no H or P division at a zero",
        },
        "all_point_component_theorem": {
            "statement": "H_r(T0)<=0 implies H_r(T)<=0 at every later point of the same connected regular physical-time interval; strict negativity remains strict",
            "selection": "At each T* select C* through edges with P(T*)!=0; all internal locks hold on a common neighborhood; boundary nulls vanish at T*",
            "clock_kinematics": "H_i=H_r/y_i; y_i'=(c_i-y_i)H_r; D_i=c_i^-1 D_r",
            "fixed_subset_identity": "-2 A_C H_r'+H_r (A_C')_fixed_C=N_C, A_C=sum_C G_i y_i²>=G_r>0, N_C=sum_C c_i y_i³ n_i>=0",
            "bounded_comparison": "H_r'=lambda_C H_r-eta_C, eta_C=N_C/(2A_C)>=0; |lambda_C|<=max_all_i |D_r log(a_i/a_r)| locally bounded on compact regular intervals",
            "positive_part_proof": "Y=max(H_r,0) is AC with Y'<=M_*Y a.e.; Gronwall from Y(T0)=0 gives sign preservation; W=-H_r comparison also preserves strict negativity",
            "algebraic_strata": "Includes isolated/persistent/multiple/accumulated roots, identically zero P and arbitrary branch histories without differentiating membership or requiring finite switches",
            "not_a_global_all_vertex_K": "Only the fixed selected subset is differentiated; actual mixed-algebraic control rejects the full-vertex kinetic sum",
            "optional_integrating_factor": "Bounded measurable lambda_C(T) gives exp(-integral lambda dT) H_r nonincreasing; a diagnostic, not a physical or canonical field map",
        },
        "actual_controls": {
            "nonempty_expanding_family": "Any finite tree with beta=(-3,1,0,0,-1), aligned a=T^(1/3),N=1,separate free phi_i=sqrt(2G_i/3)logT solves all equations; checked G=(1,0,2,3) branched orientation",
            "mixed_algebraic_family": "G=(1,1,1),a_i=T^(1/3),N=(1,1/(3T),1),edges1->0 beta=(0,1,-1/2,0,1/2),2->0 tuned beta1; separate phi0,phi2=sqrt(2/3)logT,phi1 constant",
            "actual_selected_subset": "P=(0,2),C_root=(0,2),A_C=2,N_C=4/(3T²); erroneous A_all=3 leaves residual +2/(3T²)",
            "zero_target_outside_domain": "G=(0,1),all beta/source zero,a=(1+T²,1),N=1 solves all equations with undetermined target H'(0)=2",
            "NEC_outside_domain": "Aligned a=1+T² with tuned edges and separate conserved fluids rho=3H²,p=-2H'-rho has n(0)=-4; not positive canonical matter",
        },
        "conditional_actual_CD_matching": {
            "necessary_budget": "With the actual vertex-metric/proper-time/endpoint dictionary, max_endpoint |H_r-H_CD| >=4L/[tau(1+L²)] at T=+-L*tau",
            "half_window": "L=1/2 gives 8/(5tau); errors strictly below this at both endpoints force a forbidden sign change",
            "sharpness": "Equality remains inconclusive and is attained by a tuned-beta1 static Minkowski tree with zero sources",
            "not_supplied": "No original DHOST/free-matter operator map, composite-frame identification, approximate-action-to-parent solution theorem or physical vacuum-to-bounce trajectory is supplied or assumed",
        },
        "primary_sources": {"reciprocal_edges_and_tree_leaf_stripping_prior_art": "https://arxiv.org/pdf/2410.10976",
                            "cosmological_edge_balance_and_branch_dictionary": "https://arxiv.org/pdf/2304.09205"},
        "verification_boundary": "Exact symbolic and independent Fraction/coefficient/first-jet/Gaussian replay plus separately authored covariant audits support the written all-point real-analysis proof; no proof-assistant formalization or numerical enumeration of branch histories",
        "not_established": [
            "No theorem for cycles, infinite graphs, nonpairwise or derivative interactions, or variable coefficients",
            "No theorem for shared/nonconserved/non-NEC matter, negative Einstein coefficients or a zero-EH physical target without an additional positive-component guarantee",
            "No continuation through singular/nonpositive root/lapse domains and no nonflat or anisotropic extension",
            "No scalar/vector/tensor stability, ghost, quantum, UV, positivity or cutoff verdict",
            "No universal parent exclusion, original C/D operator matching, or closure of S6, P8(b) or P8",
        ],
    }


def validate_report(expected, actual):
    if expected != actual:
        raise ValueError("Finite HR-tree sign certificate differs from exact read-only replay")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    actual = build_report()
    if args.check:
        validate_report(json.loads(REPORT.read_text()), actual)
        print("P8 S6.17.TREE: finite HR-tree physical sign-preservation replay passed; original P8 OPEN")
    else:
        print(json.dumps(actual, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
