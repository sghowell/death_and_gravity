"""Full six-invariant covariant Hessians and corrected original finite part."""

from functools import cache

import sympy as s
from p8_vacuum_affine_dimensional_spatial_symbol import matching as old_matching

from . import density, geometry, jets


def hessians(T, V, W, S=0, UD=0, UG=0):
    T, V, W, S, UD, UG = map(s.sympify, (T, V, W, S, UD, UG))
    d, t, a, H, p = geometry.d, jets.t, jets.a, jets.H, jets.p
    Hp = s.diff(H, t)
    g0, g1, g2 = s.symbols("Gamma0 Gamma1 Gamma2", real=True)
    AF, PF = a ** (d - 2), a ** (d - 4)
    E = -T + 2 * V + S - UD - UG
    AD, AG = S - UD, S - UG
    R = AF * p**2 * E * g0 / 2
    Rbackground = d * (2 * Hp + (d + 1) * H**2)
    R2 = (
        2 * Rbackground * R
        + 2
        * AF
        * p**2
        * (
            (AD + AG) * g2
            + ((d - 5) * AD + (d + 1) * AG) * H * g1
            - 3 * AD * (Hp + (d - 2) * H**2) * g0
        )
        + 2 * PF * p**4 * (S - UD - UG + W) * g0
    )
    A, B, C = T - 2 * V + UG, T - 2 * V + UD, V - UD - UG + S
    Ric = (
        AF
        * p**2
        * (
            (E * (Hp + d * H**2) - (A + AD) * (Hp + (d - 2) * H**2)) * g0
            + ((d - 4) * A / 2 + d * B / 2 + (d - 2) * C - AD + AG) * H * g1
            + ((A + B) / 2 + C) * g2
        )
        + PF * p**4 * (T - 2 * V + 2 * W + S - UD - UG) * g0 / 2
    )
    WA = -2 * A / (d - 1) + 4 * AD / (d * (d - 1))
    WC = -2 * B / (d - 1) + 4 * AG / (d * (d - 1))
    WB = -4 * (T - V) + 4 * C / (d - 1)
    W4 = (
        2 * (d - 2) * (T - 2 * V) / (d - 1)
        + 2 * (d - 2) * W / d
        - 2 * (d - 2) * (S - UD - UG) / (d * (d - 1))
    )
    Weyl = (
        AF
        * p**2
        * (
            (WA - WB + WC) * g2
            + ((2 * d - 5) * WA - (d - 2) * WB + WC) * H * g1
            + (d - 3) * WA * (Hp + (d - 2) * H**2) * g0
        )
        + PF * p**4 * W4 * g0
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
    invariants = s.symbols("T V W S UD UG", real=True)
    g, rows = hessians(*invariants)
    pole = (
        jets.m**2 * rows["R_old"]
        - rows["R_squared"] / 4
        + s.Rational(29, 30) * rows["Ricci_squared"]
        - s.Rational(2, 15) * rows["Riemann_squared"]
    )
    return invariants, g, s.factor(pole)


@cache
def invariant_logarithmic_jets():
    rows = density.invariant_spatial_coefficients()
    return tuple(
        tuple(
            (
                s.factor(
                    sum(rows[j, r, 4 - j][i] for j in range(r, 5)).subs(geometry.d, 3)
                ),
                s.factor(
                    s.diff(
                        sum(rows[j, r, 4 - j][i] for j in range(r, 5)), geometry.d
                    ).subs(geometry.d, 3)
                ),
            )
            for i in range(6)
        )
        for r in range(5)
    )


@cache
def finite_input():
    inv, g, pole = fixed_pole()
    logs = invariant_logarithmic_jets()
    zero = [
        sum(pair[0] * I for pair, I in zip(row, inv)) / (2 * s.pi**2)
        for row in logs[:3]
    ]
    slope = [
        sum(pair[1] * I for pair, I in zip(row, inv)) / (2 * s.pi**2)
        for row in logs[:3]
    ]
    power = density.invariant_spatial_coefficients()[0, 0, 2]
    F2 = sum(v.subs(geometry.d, 3) * I for v, I in zip(power, inv)) / (2 * s.pi**2)
    ell = s.Symbol("ell", real=True)
    actual = (
        sum(((1 - s.log(2) - ell / 2) * zero[r] - slope[r]) * g[r] for r in range(3))
        - jets.m**2 * F2 * g[0] / 2
        + s.diff(pole, geometry.d).subs(geometry.d, 3) / (32 * s.pi**2)
    )
    return inv, g, ell, tuple(s.factor(s.diff(actual, v)) for v in g)


@cache
def preliminary_checks():
    inv, g, pole = fixed_pole()
    logs = invariant_logarithmic_jets()
    checks = {
        "full_six_invariant_physical_covariant_pole": s.factor(
            pole.subs(geometry.d, 3)
            - 16
            * sum(
                pair[0] * I * g[r]
                for r, row in enumerate(logs[:3])
                for pair, I in zip(row, inv)
            )
        )
    }
    _, new_hess = hessians(*inv[:3])
    _, old_hess = old_matching.hessians(*inv[:3])
    for key in old_hess:
        checks["tracefree_covariant_Hessian_" + key] = s.factor(
            new_hess[key] - old_hess[key]
        )
    return checks


@cache
def finite_coefficient_norm():
    inv, _g, ell, rows = finite_input()
    L, B = s.symbols("L B", real=True)
    out, checks = {}, {}
    for r, row in enumerate(rows):
        for i, I in enumerate(inv):
            value = s.diff(row, I).subs({ell: 0, jets.m: 1000})
            polynomial = s.cancel(
                s.expand_log(value * s.pi**2 * jets.a, force=True).subs(
                    {s.log(jets.t**2 + 1): L, s.log(2): B}
                )
            )
            poly = s.Poly(polynomial, jets.t, jets.p, L, B)
            if (
                any(c.is_Rational is not True for c in poly.coeffs())
                or poly.degree(jets.p) > 4
            ):
                raise ValueError(
                    "Require the exact degree-four fixed-prescription log polynomial"
                )
            weight = (1, 1, 1, 3, 2, 2)[i]
            bound = (
                weight
                * sum(
                    abs(c)
                    * s.Rational(1, 2) ** powers[0]
                    * s.Rational(1, 4) ** powers[2]
                    for powers, c in poly.terms()
                )
                / 9
            )
            out[f"source_{r}_{I}"] = bound
            checks[f"source_{r}_{I}_full_log_polynomial_reconstruction"] = s.factor(
                polynomial - poly.as_expr()
            )
    return {"rows": out, "sum": sum(out.values()), "checks": checks}


@cache
def data():
    checks = dict(preliminary_checks())
    inv, _g, _ell, rows = finite_input()
    swap = {inv[4]: inv[5], inv[5]: inv[4]}
    transpose = [v.xreplace(swap) for v in rows]
    green = (
        rows[2] - transpose[2],
        rows[1] - 2 * s.diff(transpose[2], jets.t) + transpose[1],
        rows[0]
        - transpose[0]
        + s.diff(transpose[1], jets.t)
        - s.diff(transpose[2], jets.t, 2),
    )
    for r, v in enumerate(green):
        checks[f"full_ordered_finite_Green_identity_{r}"] = s.factor(v)
    old_rows = old_matching.finite_input()[3]
    correction = (
        jets.p**2 * (3 * jets.t**2 + 1) * (-31 * inv[0] + 30 * inv[1]) / (420 * s.pi**2)
    )
    for r, row in enumerate(rows):
        checks[f"explicit_tracefree_finite_erratum_source_jet_{r}"] = s.factor(
            row.subs(dict.fromkeys(inv[3:], 0))
            - old_rows[r]
            - (correction if r == 0 else 0)
        )
        checks[f"finite_UV_difference_zero_transfer_{r}"] = s.factor(
            row.subs(jets.p, 0)
        )
    _, curv = hessians(*inv)
    checks["complete_physical_Euler_total_derivative"] = s.factor(
        curv["Euler"].subs(geometry.d, 3)
    )
    norm = finite_coefficient_norm()
    checks.update(norm["checks"])
    checks["full_exact_Frobenius_weighted_coefficient_bound"] = norm[
        "sum"
    ] - s.Rational(66786194445073, 464486400)
    return {
        "full_covariant_Hessians": "Derive R,R^2,Ricci^2 and Riemann^2 for arbitrary trace, including the full exponential volume. EH spatial invariant is a^(d-2)p^2[-T+2V+S-UD-UG]/2. The Ricci derivation uses the conformal ADM matrix and Weyl uses the linear flat curvature with its a^(d-3) weight. All ordered time integrations precede the physical clock conversion.",
        "fixed_pole": "Keep exactly m^2 H_R-H_R2/4+29 H_Ric2/30-2 H_Riem2/15. Its full physical six-invariant value equals16 times the corrected normalized logarithmic mode density. No new finite gravitational coefficient is chosen.",
        "original_finite_part": "For fixed physical invariants, Ffinite=-m^2 F2 Gamma/2+(1-log2-ell/2)F4 Gamma-(partial_d f4 at3)/(2pi^2)+(partial_d H_pole at3)/(32pi^2). The same original radial MSbar normalization, fixed comoving splitm, Euler continuation and full volume logarithm are retained.",
        "finite_tracefree_erratum": correction,
        "ordered_finite_coefficients": rows,
        "coefficient_norm": norm,
        "norm_scope": "At m1000,ell0,|t|<=1/2 use |T|,|V|,|W|<=Frobenius product, |S|<=3 product and |UD|,|UG|<=sqrt3 product<2 product. Exact weighted log-polynomial sum66786194445073/464486400<1e6 bounds this local spatial UV difference by1e6 ||D||L2 Z24[Gamma]. It is not the complete scalar or corrected tracefree current norm.",
        "checks": checks,
        "gates": {
            "all_six_physical_covariant_pole_invariants": True,
            "same_original_finite_radial_normalization": True,
            "continued_volume_and_Euler_not_dropped": True,
            "ordered_formal_Green_not_symmetric_retarded_kernel": True,
            "nonzero_tracefree_evanescent_phase_correction": correction != 0,
            "complete_local_spatial_UV_difference_bound": norm["sum"] < 10**6,
            "full_scalar_state_time_contact_anchor_not_claimed": True,
            "S213_current_and_S215_known_input_still_withdrawn": True,
        },
    }
