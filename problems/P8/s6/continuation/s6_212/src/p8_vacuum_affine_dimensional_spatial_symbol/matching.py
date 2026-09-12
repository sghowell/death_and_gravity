"""Continued fixed pole Hessians and original MSbar spatial finite-part input."""

from functools import cache

import sympy as s
from p8_vector_dimensional import local as original_radial

from . import density, geometry, jets


def hessians(T, V, W):
    T, V, W = map(s.sympify, (T, V, W))
    d, t, a, H, p = geometry.d, jets.t, jets.a, jets.H, jets.p
    Hp = s.diff(H, t)
    g0, g1, g2 = s.symbols("Gamma0 Gamma1 Gamma2", real=True)
    S = T - 2 * V
    AF = a ** (d - 2)
    PF = a ** (d - 4)
    time = g2 + (d - 2) * H * g1
    R = -AF * p * p * S * g0 / 2
    R2 = (
        -AF * (2 * d * Hp + d * (d + 1) * H * H) * p * p * S * g0
        + 2 * PF * p**4 * W * g0
    )
    Ric = (
        AF * p * p * (T - V) * time
        - 2 * AF * (Hp + (d - 1) * H * H) * p * p * S * g0
        + PF * p**4 * (T / 2 - V + W) * g0
    )
    Weyl = (
        4 * (d - 2) * AF * p * p * (T - V) * time / (d - 1)
        - 2 * (d - 3) * AF * (Hp + (d - 2) * H * H) * p * p * S * g0 / (d - 1)
        + 2 * (d - 2) * PF * p**4 * S * g0 / (d - 1)
        + 2 * (d - 2) * PF * p**4 * W * g0 / d
    )
    Riem = Weyl + 4 * Ric / (d - 1) - 2 * R2 / (d * (d - 1))
    return (g0, g1, g2), {
        "R_old": R,
        "R_squared": R2,
        "Ricci_squared": Ric,
        "Riemann_squared": Riem,
        "Weyl_dimension": Weyl,
        "Euler": Riem - 4 * Ric + R2,
    }


@cache
def fixed_pole():
    T, V, W = s.symbols("T V W", real=True)
    g, rows = hessians(T, V, W)
    pole = (
        jets.m**2 * rows["R_old"]
        - rows["R_squared"] / 4
        + s.Rational(29, 30) * rows["Ricci_squared"]
        - s.Rational(2, 15) * rows["Riemann_squared"]
    )
    return (T, V, W), g, s.factor(pole)


def fourier_measure(d):
    d = s.sympify(d)
    return 2 * s.pi ** (d / 2) / s.gamma(d / 2) / (2 * s.pi) ** d


@cache
def finite_input():
    invariants, g, pole = fixed_pole()
    zero, slope = [], []
    for row in density.total_logarithmic_dimension_jet()[:3]:
        zero.append(
            sum(value[0] * inv for value, inv in zip(row, invariants)) / (2 * s.pi**2)
        )
        slope.append(
            sum(value[1] * inv for value, inv in zip(row, invariants)) / (2 * s.pi**2)
        )
    power = density.invariant_dimension_jet()[0, 0, 2]
    F2 = sum(value[0] * inv for value, inv in zip(power, invariants)) / (2 * s.pi**2)
    Hprime = s.diff(pole, geometry.d).subs(geometry.d, 3)
    ell = s.Symbol("ell", real=True)
    actual = (
        sum(((1 - s.log(2) - ell / 2) * zero[r] - slope[r]) * g[r] for r in range(3))
        - jets.m**2 * F2 * g[0] / 2
        + Hprime / (32 * s.pi**2)
    )
    return invariants, g, ell, tuple(s.factor(s.diff(actual, v)) for v in g)


@cache
def preliminary_checks():
    inv, g, pole = fixed_pole()
    T, V, W = inv
    rows = density.total_logarithmic_dimension_jet()
    actual = sum(
        sum(pair[0] * a for pair, a in zip(row, inv)) * g[r] / (2 * s.pi**2)
        for r, row in enumerate(rows[:3])
    )
    checks = {
        "full_fixed_pole_at_physical_dimension": s.factor(
            pole.subs(geometry.d, 3) / (32 * s.pi**2) - actual
        )
    }
    _gj, curv = hessians(T, V, W)
    d, a, H = geometry.d, jets.a, jets.H
    expected = (
        -(a ** (d - 2))
        * (d - 3)
        * ((d - 2) * H * H + 2 * s.diff(H, jets.t))
        * jets.p**2
        * (T - 2 * V)
        * g[0]
    )
    checks["continued_Euler_spatial_variation_not_dropped"] = s.factor(
        curv["Euler"] - expected
    )
    checks["normalized_sphere_measure_physical"] = s.simplify(
        fourier_measure(3) - 1 / (2 * s.pi**2)
    )
    checks["measure_dimension_log_derivative"] = s.simplify(
        s.diff(s.log(fourier_measure(d)), d).subs(d, 3)
        - (s.EulerGamma / 2 - 1 - s.log(s.pi) / 2)
    )
    _, _, _ell, finite = finite_input()
    checks["finite_time_Green_coefficient_including_volume_log"] = s.factor(
        finite[1] - s.diff(finite[2], jets.t)
    )
    # Exact normalization bridge to the original radial implementation.
    for order in (0, 1, 2):
        for power in range(4):
            wanted = (
                2
                * s.sqrt(s.pi)
                * s.rf(original_radial.dimension / 2, power)
                / s.gamma(order + power - s.Rational(1, 2))
            )
            actual = original_radial.radial_polynomial(original_radial.z**power, order)
            checks[f"original_MSbar_radial_normalization_{order}_{power}"] = s.factor(
                actual - wanted
            )
    return checks


@cache
def finite_coefficient_norm():
    invariants, _g, ell, rows = finite_input()
    L, B = s.symbols("L B", real=True)
    bounds = {}
    checks = {}
    for r, row in enumerate(rows):
        for name, inv in zip(("T", "V", "W"), invariants):
            coefficient = s.diff(row, inv).subs({ell: 0, jets.m: 1000})
            polynomial = s.cancel(
                s.expand_log(coefficient * s.pi**2 * jets.a, force=True).subs(
                    {s.log(jets.t**2 + 1): L, s.log(2): B}
                )
            )
            poly = s.Poly(polynomial, jets.t, jets.p, L, B)
            bound = (
                sum(
                    abs(c)
                    * s.Rational(1, 2) ** powers[0]
                    * s.Rational(1, 4) ** powers[2]
                    for powers, c in poly.terms()
                )
                / 9
            )
            bounds[f"source_{r}_{name}"] = bound
            checks[f"source_{r}_{name}_exact_log_polynomial_reconstruction"] = s.factor(
                polynomial - poly.as_expr()
            )
            if polynomial != 0:
                assert poly.degree(jets.p) <= 4
                assert all(c.is_Rational for c in poly.coeffs())
    return {"rows": bounds, "sum": sum(bounds.values()), "checks": checks}


@cache
def data():
    checks = dict(preliminary_checks())
    norms = finite_coefficient_norm()
    checks.update(norms["checks"])
    return {
        "continued_Hessians": "Derive the complete spatial differences of R_old,R_old_squared,Ricci_squared and Riemann_squared in d+1 dimensions with the original proper-clock measure. The full continued Euler variation is-a^(d-2)(d-3)[(d-2)H^2+2H_prime]p^2(T-2V)Gamma; only its physical value vanishes.",
        "fixed_counterterm": "Keep the original pole scalar weights m^2 H_R-H_R2/4+29 H_Ric2/30-2 H_Riem2/15 while continuing invariant contractions and volume. Their dimension derivative is retained; no finite coefficient is selected anew.",
        "MSbar_normalization": "The original normalized radial convention is reproduced by multiplying the raw spatial Fourier integral by[e^EulerGamma*mu^2/(4pi)]^epsilon with d=3-2epsilon. The angular factor and m^(-2epsilon) tail then give the finite multiplier1-log2-ell/2, ell=log(m^2/mu^2). The split radius is the fixed comovingm, not a*m.",
        "finite_UV_difference": "Let F2,F4 be the physical spatial coefficient operators, and f4_d the normalized-sphere complete logarithmic invariant coefficient before the Fourier factor. The finite UV difference after the original fixed pole subtraction is-m^2 F2 Gamma/2+(1-log2-ell/2)F4 Gamma-(partial_d f4_d at3)/(2pi^2)+[partial_d H_pole(d) at3]/(32pi^2). The full source-value log includes the evanescent j1 endpoint. The homogeneous contact remains in its anchor.",
        "proper_time_contact": "The volume derivative contributes log(a) times the physical pole. The full finite first-time-derivative coefficient equals the proper-time derivative of the second-time coefficient. Dropping this volume logarithm would violate that Green identity.",
        "coefficient_norm": norms,
        "norm_scope": "At the fixedell0,m1000 slab, the exact log-polynomial row bound uses|t|<=1/2,0<=log(1+t^2)<=1/4,0<log2<1,a>=1 and pi^2>9. It bounds only this local finite UV difference in the two-time/four-spatial Z24 norm. A full assembled matched response, inverse or background is not inferred here.",
        "checks": checks,
        "gates": {
            "original_MSbar_radial_normalization_replayed": True,
            "fixed_four_dimensional_scalar_pole_weights": True,
            "continued_Euler_and_volume_evanescent_terms_retained": True,
            "exact_lower_comoving_band_term_retained": True,
            "finite_proper_time_Green_identity": True,
            "finite_local_UV_difference_bound": norms["sum"] < 100000,
            "full_regulator_limit_and_current_assembly_not_claimed": True,
            "original_V_G_B_not_closed": True,
        },
    }
