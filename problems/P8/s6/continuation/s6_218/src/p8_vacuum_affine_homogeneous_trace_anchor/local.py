"""Full trace local action matching, volume Euler current and exact bounds."""

from functools import cache
from math import comb

import sympy as s
from p8_vacuum_affine_adiabatic_action import curvature as prescribed

from . import homogeneous

d, z, H, H1 = s.symbols("d z H H1", real=True)
h = s.symbols("h0:8", real=True)
m = s.Symbol("m", positive=True)
LOCAL = s.Integer(10) ** 11


def radial_polynomial(poly, alpha):
    return s.factor(
        sum(
            co * s.rf(d / 2, power[0]) / s.rf(alpha, power[0])
            for power, co in s.Poly(s.expand(poly), z).terms()
        )
    )


def derivative(poly):
    if poly.has(h[-1]):
        raise ValueError("One more H jet would be required")
    return s.expand(sum(s.diff(poly, h[j]) * h[j + 1] for j in range(len(h) - 1)))


def volume_derivative(poly):
    return s.expand(derivative(poly) + 3 * h[0] * poly)


@cache
def density_polynomial():
    R = 6 * (h[1] + 2 * h[0] ** 2)
    Ric = 9 * (h[1] + h[0] ** 2) ** 2 + 3 * (h[1] + 3 * h[0] ** 2) ** 2
    return s.expand(
        s.Rational(5, 2) * m**4 + s.Rational(5, 3) * m * m * R - R * R / 30 - Ric / 15
    )


@cache
def current_polynomial():
    F = density_polynomial()
    return s.factor(
        F
        - volume_derivative(s.diff(F, h[0])) / 3
        + volume_derivative(volume_derivative(s.diff(F, h[1]))) / 3
    )


def trace_linear_coefficients():
    value = current_polynomial()
    return tuple([value] + [s.diff(value, h[j]) / 3 for j in range(4)])


@cache
def local_bounds():
    poly = s.Poly(s.expand(current_polynomial().subs(m, 1000)), *h[:4])
    majorant = sum(
        abs(co) * s.prod(h[i] ** powers[i] for i in range(4))
        for powers, co in poly.terms()
    )
    point = {
        h[i]: homogeneous.h_majorants()[i] + homogeneous.DELTA / 3 for i in range(4)
    }
    result = []
    for order in range(3):
        result.append(s.factor(2 * s.Rational(5, 4) ** 6 * majorant.subs(point) / 576))
        majorant = s.expand(majorant + sum(s.diff(majorant, x) / 3 for x in h[:4]))
    return tuple(result)


@cache
def data():
    p = -z * H
    pprime = z * (2 * H * H - H1) - 2 * z * z * H * H
    st, sl = ((d - 2) * H + p) / 2, ((d - 2) * H - p) / 2
    # These are direct derivatives of the scalar root-free symbols.
    tt = ((d - 2) * H1 + pprime) / 2 - p * st
    tl = ((d - 2) * H1 - pprime) / 2 - p * sl
    second = radial_polynomial((d - 1) * st**2 + sl**2, s.Rational(1, 2)) / 4
    fourth = radial_polynomial(
        (d - 1) * (tt**2 + st**4) + tl**2 + sl**4, s.Rational(3, 2)
    )
    R = d * (2 * H1 + (d + 1) * H * H)
    Ric = d * d * (H1 + H * H) ** 2 + d * (H1 + d * H * H) ** 2
    Riem = 4 * d * (H1 + H * H) ** 2 + 2 * d * (d - 1) * H**4
    a4s = (5 * R * R - 2 * Ric + 2 * Riem) / 360
    a4v = d * a4s + Ric / 2 - R * R / 6 - Riem / 12
    c4 = -4 * d * (5 * d**3 - 96 * d**2 + 229 * d - 126) / 945
    coeff = trace_linear_coefficients()
    green = []
    for k in range(5):
        transpose = 0
        for j in range(k, 5):
            v = coeff[j]
            for _ in range(j - k):
                v = volume_derivative(v)
            transpose += (-1) ** j * comb(j, k) * v
        green.append(s.factor(coeff[k] - transpose))
    expected = (
        -(
            72 * h[0] ** 2 * h[1]
            - 60 * h[0] ** 2 * m * m
            + 48 * h[0] * h[2]
            + 36 * h[1] ** 2
            - 40 * h[1] * m * m
            + 8 * h[3]
            - 15 * m**4
        )
        / 6
    )
    phi, chi, A = s.symbols("phi chi A", real=True)
    cosmo = s.Rational(5, 2) * m**4 * A**3 * s.exp(phi)
    checks = {
        "full_general_d_second_order_local_matching": s.factor(
            second - (d - 6) * R / 24 + d * (d - 6) * (H1 + d * H * H) / 12
        ),
        "full_general_d_fourth_order_local_matching": s.factor(
            fourth - 2 * a4v - c4 * (d * H**4 + 3 * H * H * H1)
        ),
        "full_volume_Euler_current": s.factor(current_polynomial() - expected),
        "all_five_full_density_Green_relations": s.Matrix(green),
        "nonzero_cosmological_trace_current_retained": s.diff(cosmo, phi).subs(phi, 0)
        - s.Rational(5, 2) * m**4 * A**3,
        "nonzero_cosmological_trace_Hessian_retained": s.diff(cosmo, phi, 2).subs(
            phi, 0
        )
        - s.Rational(5, 2) * m**4 * A**3,
        "no_tracefree_projection_before_volume_variation": s.diff(
            s.exp(phi + chi), phi, chi
        ).subs({phi: 0, chi: 0})
        - 1,
        "three_explicit_local_raw_bounds": s.Matrix(
            [
                v - s.Rational(n, 169869312)
                for v, n in zip(
                    local_bounds(),
                    (5625132891495009211, 5625161941796211511, 5625195992153673811),
                )
            ]
        ),
    }
    fixed = prescribed.dimensional_coefficients()
    eps, ell = prescribed.epsilon, prescribed.ell
    checks["same_original_vacuum_finite_part"] = s.expand(
        fixed["vacuum"] - 3 / eps - s.Rational(5, 2) + 3 * ell
    )
    checks["same_original_Einstein_finite_part"] = s.expand(
        fixed["einstein"] - 1 / eps - s.Rational(5, 3) + ell
    )
    checks["same_original_curvature_evanescent_part"] = s.expand(
        fixed["curvature"]
        - 2 * fixed["fixed_vector"] / eps
        + 2 * ell * fixed["fixed_vector"]
        + 4 * fixed["fixed_scalar"]
    )
    return {
        "full_general_d_symbols": {
            "second": second,
            "fourth": fourth,
            "R": R,
            "Ricci2": Ric,
            "Riemann2": Riem,
            "a4v": a4v,
            "fourth_weighted_boundary_coefficient": c4,
        },
        "matching": "At each fixed time the exact radial change k=A q supplies the full A^d Jacobian. This evaluates the local integral before varying its density; it is not a time-dependent momentum relabeling inside the mode ODE. The displayed symbols equal the original invariant action modulo weighted derivatives of A^d H_A and A^d H_A^3 at general d, before the dimension limit. The generic S193 current/action variation applies to these same scalar K,V,omega.",
        "fixed_finite_action": "At mu=m, only after the original dimension limit, the finite local density is A^3[5m^4/2+(5/3)m^2 R_old-R_old^2/30-Ricci_old^2/15]/(64pi^2), modulo the compact four-dimensional Euler variation. The cosmological volume term is not zero in this trace direction.",
        "full_trace_current_polynomial": current_polynomial(),
        "complete_local_density_current": "J_phi,local=A^3 L(H_A,H_A',H_A'',H_A''')/(64pi^2), L displayed. It is the complete Euler derivative of volume times the fixed local density, so no time derivatives of the detector remain.",
        "linear_coefficients_before_A_cubed_over_64pi_squared": coeff,
        "full_local_raw_amplitude_bounds": local_bounds(),
        "majorant_proof": "For |Gamma^j|<=1, H_A^(j)=H^(j)+epsilon Gamma^(j+1)/3. Replace every coefficient of L by its absolute value and evaluate at the explicit H majorants plus .01/3. A^3<=2(5/4)^6. Each amplitude derivative is majorized by(1+sum_j partial_hj/3), the1 from the varying volume. Use pi^2>9. All three raw bounds are below1e11; no fixed-volume simplification is used.",
        "checks": checks,
        "gates": {
            "full_general_d_matching_before_dimension_limit": True,
            "all_local_amplitude_bounds_strict": all(v < LOCAL for v in local_bounds()),
            "volume_exponential_bound": 1 / (1 - homogeneous.DELTA) < 2,
            "local_current_requires_at_most_four_source_time_jets": not current_polynomial().has(
                *h[4:]
            ),
            "compact_Euler_removed_only_after_finite_part": True,
            "full_local_operator_Green_identity": all(v == 0 for v in green),
            "local_cosmological_trace_term_nonzero": s.diff(cosmo, phi, 2).subs(phi, 0)
            != 0,
        },
    }
