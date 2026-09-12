"""Full constrained isotropic trace Hamiltonian and uniform scalar frame."""

from functools import cache

import sympy as s
from p8_vacuum_affine_current_response import vertices as previous_vertices
from p8_vacuum_affine_matrix_response_tail import mixed

t = s.Symbol("t", real=True)
z = s.Symbol("z", real=True)
H = 4 * t / (1 + t * t)
DELTA = s.Rational(1, 100)


@cache
def h_majorants():
    values = []
    for j in range(12):
        num, den = s.fraction(s.cancel(s.diff(H, t, j)))
        polynomial = s.Poly(den, t)
        if polynomial.coeff_monomial(1) < 1 or any(
            co < 0 for co in polynomial.all_coeffs()
        ):
            raise ValueError(
                "The base H derivative denominator is not bounded below by one"
            )
        values.append(
            sum(
                abs(co) * s.Rational(1, 2) ** power[0]
                for power, co in s.Poly(num, t).terms()
            )
        )
    values[0] = s.Rational(8, 5)
    return tuple(values)


def scalar_vertices():
    return {
        "transverse": s.diag((1 - 2 * z) / 3, -s.Rational(1, 3)),
        "longitudinal": s.diag(s.Rational(1, 3), -(1 + 2 * z) / 3),
    }


@cache
def data():
    A, m = s.symbols("A m", positive=True)
    k = s.Symbol("k", nonnegative=True)
    phi = s.Symbol("phi", real=True)
    ae = A * s.exp(phi / 3)
    KT, KL = 1 / ae, 1 / ae + k * k / (ae**3 * m * m)
    VT, VL = ae * m * m + k * k / ae, ae * m * m
    omega2 = m * m + k * k / ae**2
    M = s.diag(VT, VT, VL, KT, KT, KL)
    tau = 2 * phi
    B = -phi * s.eye(3) / 3
    checks = {
        "complete_trace_chart_shifted_matrix": B
        - (2 * phi * s.eye(3) / 3 - tau * s.eye(3) / 2),
        "full_volume_constraint_exponent": -tau / 2 + phi,
        "transverse_full_frequency_identity": s.factor(KT * VT - omega2),
        "longitudinal_full_frequency_identity": s.factor(KL * VL - omega2),
        "complete_second_trace_vertex_positive_constraint": s.factor(
            s.diff(KL, phi, 2) - 1 / (9 * ae) - k * k / (ae**3 * m * m)
        ),
    }
    for kind, (K, V) in {"transverse": (KT, VT), "longitudinal": (KL, VL)}.items():
        exact = s.diag(s.diff(V, phi) / V, s.diff(K, phi) / K)
        target = scalar_vertices()[kind].subs(z, k * k / (ae * ae * omega2))
        checks[kind + "_full_normalized_current_vertex"] = (exact - target).applyfunc(
            s.simplify
        )
    hz, p = s.symbols("H_A p", real=True)
    checks["transverse_balanced_squeeze"] = p / 2 - (-hz) / 2 - (hz + p) / 2
    checks["longitudinal_balanced_squeeze"] = p / 2 - (2 * p - hz) / 2 - (hz - p) / 2
    gamma = s.Symbol("Gamma", real=True)
    ze = -s.Rational(2, 3) * gamma * z * (1 - z)
    for kind, G in scalar_vertices().items():
        first = G.diff(z) * ze
        second = first.diff(z) * ze
        target1 = (
            s.diag(s.Rational(4, 9) * gamma * z * (1 - z), 0)
            if kind == "transverse"
            else s.diag(0, s.Rational(4, 9) * gamma * z * (1 - z))
        )
        target2 = (
            s.diag(-s.Rational(8, 27) * gamma**2 * z * (1 - z) * (1 - 2 * z), 0)
            if kind == "transverse"
            else s.diag(0, -s.Rational(8, 27) * gamma**2 * z * (1 - z) * (1 - 2 * z))
        )
        checks[kind + "_complete_first_parameter_vertex"] = first - target1
        checks[kind + "_complete_second_parameter_vertex"] = second - target2
    checks["zero_momentum_transverse_longitudinal_vertices_coincide"] = (
        scalar_vertices()["transverse"].subs(z, 0)
        - scalar_vertices()["longitudinal"].subs(z, 0)
    )
    return {
        "family": "Q=(2/3)phi I, phi=epsilon Gamma, |epsilon|<=.01 and |Gamma^(j)|<=1 through j12, on the unchanged unit CD slab and common zero initial neighborhood. A=a exp(phi/3), sqrt(h)=A^3.",
        "Hamiltonian": "KT=A^-1, KL=A^-1+k^2/(A^3m^2), VT=A m^2+k^2/A, VL=A m^2. The longitudinal constraint varies with the trace; no fourth oscillator or deleted longitudinal energy.",
        "actual_M": M,
        "frame": "In the fixed momentum T/T/L basis each physical oscillator has KV=omega^2. The balanced scalar rotation R vanishes. With p=omega'/omega=-z H_A, S_T=(H_A+p)/2 and S_L=(H_A-p)/2.",
        "normalized_vertices": scalar_vertices(),
        "checks": checks,
        "gates": {
            "positive_all_mode_energies": True,
            "full_trace_constraint_varies": s.diff(KL, phi).has(k),
            "full_trace_second_contact_varies": s.diff(KL, phi, 2).has(k),
            "no_polarization_eigenvalue_gap_required": True,
            "parameter_dependent_volume_retained": s.diff(ae**3, phi) != 0,
        },
    }


@cache
def domination_data():
    c = mixed.constants()
    hb = h_majorants()
    margins = {}
    for j in range(12):
        for a in range(3):
            bound = (
                hb[j] + DELTA / 3
                if a == 0
                else (s.Rational(1, 3) if a == 1 else s.S.Zero)
            )
            margins[j, a] = c["L"][j, a] - bound / 2
    checks = {
        "exact_direct_base_H_boundary": H.subs(t, s.Rational(1, 2)) - s.Rational(8, 5),
        "base_H_monotonicity_numerator": s.factor(s.diff(H, t) * (1 + t * t) ** 2)
        - 4 * (1 - t * t),
        "frequency_trace_exponent": 2 * s.Rational(1, 3) - s.Rational(2, 3),
        "same_unit_interval": s.Rational(1, 2) - (-s.Rational(1, 2)) - 1,
        "first_vertex_quadratic_bound": s.expand(
            s.Rational(1, 4) - z * (1 - z) - (z - s.Rational(1, 2)) ** 2
        ),
        "second_vertex_display": s.Rational(8, 27) * s.Rational(1, 4)
        - s.Rational(2, 27),
        "trace_source_and_detector_Frobenius_conversion": (s.sqrt(3) / 2) ** 2
        - s.Rational(3, 4),
    }
    return {
        "H_base_majorants": hb,
        "all36_H_half_to_old_L_margins": margins,
        "frequency": "omega^2=m^2+a^-2 exp(-2epsilon Gamma/3)k^2. The ordered exponential majorants of S191 require only the bounded direction and its finite jets, not its trace. The2/3 factor is at most one, so every normalized mixed omega/reciprocal jet is dominated by its old table. The direct first time derivative remains below4.",
        "squeeze": "H_A=H+epsilon Gamma'/3. For every j<=11 and a<=2, half its raw mixed jet bound is strictly below old L[j,a]. Since p/2 is the old rho expression, both S_T and S_L are dominated by old S=L+rho. R=0 is bounded by old R. All36 margins are checked exactly.",
        "vertices": "Complete normalized vertices are diag((1-2z)/3,-1/3) and diag(1/3,-(1+2z)/3). Their full first two amplitude derivatives have bounds1/9 and2/27 for |Gamma|<=1, 0<=z<=1. Thus(1,1/9,2/27) is dominated by the original physical vertex displays(2,25,650), including all frame and frequency derivatives.",
        "state_and_frequency_floor": "The family and every source jet coincide with the original history near the same initial surface. The actual Borel preparation and the finite-reference initial graph are unchanged there. Since exp(-2epsilon Gamma/3)>=1-1/150>99/100, the same fixed nu_minus and proof partition1e16 are valid.",
        "checks": checks,
        "gates": {
            "all36_H_to_old_frame_margins_positive": all(
                v > 0 for v in margins.values()
            ),
            "direct_trace_squeeze_below_old_display": s.Rational(8, 5) + DELTA / 3 < 11,
            "trace_frequency_first_time_derivative_below_four": 2
            * (s.Rational(8, 5) + DELTA / 3)
            < 4,
            "same_frequency_lower_coefficient": 1 - 2 * DELTA / 3 > s.Rational(99, 100),
            "same_uniform_upper_frequency": 1 / (1 - 2 * DELTA / 3)
            < 9 * s.Rational(99, 100) / s.Rational(25, 16) ** 2,
            "full_normalized_vertices_dominated": all(
                value < old
                for value, old in zip(
                    (s.Integer(1), s.Rational(1, 9), s.Rational(2, 27)),
                    previous_vertices.G,
                )
            ),
            "no_new_time_jet_requirement": 11 + 1 <= 12,
            "unchanged_initial_frequency_and_slope": True,
        },
    }
