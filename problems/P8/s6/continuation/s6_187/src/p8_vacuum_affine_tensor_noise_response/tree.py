"""Actual unimodular tensor action, physical source sign and canonical scaling."""

from functools import cache

import sympy as s
from p8 import tensor as original_tensor
from p8_affine_vacuum_domain import family
from p8_vacuum_affine_proca_gaussian.bridge import KAPPA
from p8_vacuum_analytic_affine_parent import source


@cache
def data():
    e, H = s.symbols("epsilon H", real=True)
    plus, cross, pdot, cdot = s.symbols(
        "gamma_plus gamma_cross plus_dot cross_dot", real=True
    )
    gamma = s.Matrix([[plus, cross, 0], [cross, -plus, 0], [0, 0, 0]])
    gd = s.Matrix([[pdot, cdot, 0], [cdot, -pdot, 0], [0, 0, 0]])
    unit = s.eye(3) + e * gamma + e * e * gamma * gamma / 2
    inv = s.eye(3) - e * gamma + e * e * gamma * gamma / 2
    derivative = 2 * H * unit + e * gd + e * e * (gd * gamma + gamma * gd) / 2
    K = (inv * derivative / 2).applyfunc(lambda x: s.series(x, e, 0, 3).removeO())
    scalar = (s.trace(K) ** 2 - s.trace(K * K)).expand()
    Rclock = family.data()["R"].subs(family.X, 1)
    old = original_tensor.derive()
    symbols = {str(x): x for x in old["GT"].free_symbols | old["FT"].free_symbols}
    GT = old["GT"].subs(symbols["B"], -Rclock / 2)
    FT = old["FT"].subs(symbols["f"], -Rclock / 2)
    h, hd, hz, gammad, gammaz, q, a = s.symbols(
        "h h_dot h_z gamma_dot gamma_z tensor_stress a", real=True
    )
    # Unit Frobenius polarization convention; gamma_tensor=2h_tensor/sqrt(kappa).
    original = KAPPA * (gammad * gammad - gammaz * gammaz / a**2) / 8
    canonical = original.subs(
        {gammad: 2 * hd / s.sqrt(KAPPA), gammaz: 2 * hz / s.sqrt(KAPPA)}
    )
    # In +---, delta S_m=(1/2) integral sqrt(-g) T_ab delta g^ab.
    # delta g^ij=gamma_ij/a^2 on the base, hence the source pairing is positive.
    coupling = q * (2 * h / s.sqrt(KAPPA)) / 2
    S = source.data()["shifted_source_P8"]
    names = {str(x): x for x in S.free_symbols}
    S_t = S.subs(
        {
            names["R"]: 1,
            names["positive_X"]: 1,
            names["Box_P8_u"]: 3 * H,
            names["H_clock"]: H,
            names["R_u"]: 0,
            names["uHu_P8"]: 0,
        }
    )
    checks = {
        "complete_R_on_unit_clock": Rclock - 1,
        "unimodular_tensor_determinant_through_quadratic": s.expand(unit.det()).coeff(
            e, 2
        ),
        "tensor_extrinsic_trace_first": s.expand(s.trace(K)).coeff(e, 1),
        "tensor_extrinsic_trace_second": s.expand(s.trace(K)).coeff(e, 2),
        "complete_extrinsic_tensor_kinetic": s.expand(
            -scalar.coeff(e, 2) / 2 - s.trace(gd * gd) / 8
        ),
        "full_source_zero_on_unimodular_tensor_history": s.simplify(S_t),
        "actual_tensor_GT_per_kappa": GT - 1,
        "actual_tensor_FT_per_kappa": FT - 1,
        "actual_tensor_spatial_R3": old["residuals"]["tensor_R3"],
        "canonical_unit_polarization_kinetic": s.simplify(
            canonical - (hd * hd - hz * hz / a**2) / 2
        ),
        "physical_canonical_tensor_stress_pairing": s.simplify(
            coupling - q * h / s.sqrt(KAPPA)
        ),
    }
    return {
        "metric_chart": "g_ij=-a^2(exp gamma)_ij, gamma symmetric transverse traceless; N=1, shift=0,u=t",
        "full_parent_restriction": "X=1,R=1,trace K=3H,uHu=0 and S=0 exactly. The fixed scalar coefficient and original homogeneous M1 have no tensor dependence in this unimodular chart; no old profile is substituted.",
        "quadratic_classical_tensor_action": "kappa/8 integral a^3 [tr gamma_t^2-a^-2 sum_i tr(partial_i gamma)^2]",
        "canonical_tensor": "h=sqrt(kappa)gamma/2 in unit-Frobenius transverse traceless polarizations",
        "tree_tensor_operator": "L=partial_t^2+3H partial_t-a^-2 Delta, H=4t/(1+t^2)",
        "leading_source_only_equation": "L h = Pi_TT T_centered/sqrt(kappa); this defines the tree-propagated vector-noise component, not the full quantum tensor equation",
        "self_adjoint_weight": "a^3 dt dx",
        "checks": checks,
        "gates": {
            "actual_tensor_kinetic_positive": GT > 0,
            "actual_tensor_speed_one": GT == FT,
            "same_fixed_kappa": KAPPA == 10**800,
        },
    }
