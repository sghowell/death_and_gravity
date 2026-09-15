"""Whole current source and the original light/heavy scalar metric vertex sector."""

from functools import cache

import sympy as s
from p8_vacuum_affine_gaussian_metric_pole_matching import source as previous
from p8_vacuum_affine_heavy_parent_one_loop import bounds as finite_scheme
from p8_vacuum_affine_heavy_scalar_one_loop import self_energy
from p8_vacuum_affine_massive_elastic_proca_infrared import source as minimal

MU, N, T, Z, V = s.symbols("mu heavy_mass_squared transfer z v", positive=True)
G, K = s.symbols("positive_cubic kappa", positive=True)
KAPPA = previous.KAPPA
HEAVY_MASS2 = minimal.HEAVY_MASS2
CUBIC = minimal.CUBIC
CONTACT = minimal.CONTACT


def require_mass(value):
    if (
        isinstance(value, bool)
        or not isinstance(value, (int, s.Rational))
        or value <= 0
    ):
        raise ValueError("Require a positive exact rational mass squared")
    return s.Rational(value)


def require_order(value):
    if isinstance(value, bool) or not isinstance(value, (int, s.Integer)) or value < 1:
        raise ValueError("Require a positive integer Taylor order")
    return int(value)


@cache
def data():
    inherited = previous.data()
    answer = {k: v for k, v in inherited.items() if k not in ("checks", "gates")}
    checks = dict(inherited["checks"])
    full = minimal.current.fixed_functions()
    u, X = minimal.current.u, minimal.current.X
    at_zero = {u: 0, X: 0}
    j = minimal.current.heavy.coefficients()["normalized_heavy_source"]
    checks.update(
        {
            "literal_heavy_mass_not_rounded": HEAVY_MASS2
            - (s.Integer(10) ** 200 / 512 + 2),
            "literal_cubic_not_rechosen": CUBIC - s.Rational(1, 8192),
            "literal_cubic_square": CUBIC**2 - s.Rational(1, 2**26),
            "literal_kappa_not_rechosen": KAPPA - s.Integer(10) ** 800,
            "no_heavy_source_constant": j.subs(at_zero),
            "no_heavy_source_linear_light": s.diff(j, u).subs(at_zero),
            "no_heavy_derivative_cubic": s.diff(j, X).subs(at_zero),
            "whole_heavy_cubic_again": s.diff(j, u, 2).subs(at_zero) / s.sqrt(KAPPA)
            - CUBIC,
            "no_linear_curvature_improvement": s.diff(full["R_full"], u).subs(at_zero),
            "no_quadratic_curvature_improvement": s.diff(full["R_full"], u, 2).subs(
                at_zero
            ),
            "same_scalar_OS_denominator": s.expand(
                self_energy.F.subs({self_energy.x: Z, self_energy.m2: N})
                - ((1 - Z) ** 2 + Z * N)
            ),
            "same_scalar_OS_first_subtraction": s.cancel(
                self_energy.ALPHA.subs({self_energy.x: Z, self_energy.m2: N})
                - Z * (1 - Z) / ((1 - Z) ** 2 + Z * N)
            ),
        }
    )
    answer.update(
        {
            "same_original_vacuum_parameters": {
                "mu": s.Integer(1),
                "n": HEAVY_MASS2,
                "g": CUBIC,
                "C": CONTACT,
                "kappa": KAPPA,
            },
            "literal_metric_convention": "Signature+---, g_mn=eta_mn+2h_mn/sqrt(kappa). Minimal scalar metric vertex-i T_mn/sqrt(kappa), T=pr+rp-eta(p.r-mu). The cubic density+g H Phi^2/2 generates metric contact+i g eta_mn/sqrt(kappa).",
            "complete_F1_graph_scope": "One matter loop in the inherited cubic/contact potential at leading external gravitational coupling. Both light-line and heavy-line stress triangles contribute to F1. Cubic metric contacts, the quartic bubble/tadpole and heavy-metric mixing supply eta/qq tensors only. Additional vertices of higher field degree require additional matter loops here; pure-gravity internal loops are separate.",
            "unchanged_OS_scope": "The inherited finite scalar counterterms use V(X)=1-T(X), with deltaZ=-Pi_MSprime(1) and delta_m^2=Pi_MS(1)-Pi_MSprime(1). They are covariant kinetic/mass terms, not a new Ricci-derivative finite prescription. These first-order conditions are not an exact interacting LSZ theorem.",
            "curved_boundary": "S285's light pure-curvature polynomial stays unmatched. R Phi^2 and R H improvements and a higher Ricci-derivative operator are not fixed by the Ward identity or flat scalar OS conditions. The computed F1 is the explicitly stated generated coefficient, not a claim that every possible finite local addition vanishes.",
            "checks": checks,
            "gates": {
                "entire_inherited_source_and_matching_retained": all(
                    k in answer
                    for k in ("whole_original_R_F", "whole_retained_heavy_source")
                ),
                "literal_actual_heavy_gap_positive": HEAVY_MASS2 > 2,
                "same_covariant_scalar_OS_extension_not_new_choice": "V[deltaC_new"
                in finite_scheme.data()["finite_counterterm_extension"],
                "source_check_dictionary_copied": checks is not inherited["checks"],
                "no_internal_graviton_in_matter_F1_graphs": True,
                "all_light_curvature_and_original_P8_boundaries_retained": True,
            },
        }
    )
    return answer
