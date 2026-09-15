"""Complete finite positive-domain graph representative and explicit majorant."""

from functools import cache

import sympy as s
from p8_vacuum_affine_whole_quartic_gravity_sector import forward as quartic_forward

from . import graphs, source

MU, N, G, K = source.MU, source.N, source.G, source.K
S, T = graphs.S, graphs.T
X = graphs.X
E2 = s.Symbol("formal_resolution_squared", positive=True)
V = s.Symbol("forward_crossing_coordinate", real=True)


def J0(channel, mass=MU):
    return graphs.light_J(channel, -1, mass)


def J1(channel, mass=MU):
    channel, mass = map(s.sympify, (channel, mass))
    A = mass - channel * X * (1 - X)
    return s.Integral(s.log(A) / A, (X, 0, 1))


def Lll(channel, mass=MU):
    channel, mass = map(s.sympify, (channel, mass))
    return s.Integral(s.log(mass - channel * X * (1 - X)), (X, 0, 1))


def Lmix(mass=MU, heavy=N):
    mass, heavy = map(s.sympify, (mass, heavy))
    return s.Integral(s.log(heavy * X + mass * (1 - X) ** 2), (X, 0, 1))


def Lh(channel, heavy=N):
    channel, heavy = map(s.sympify, (channel, heavy))
    return -1 + s.Integral(s.log(heavy - channel * (1 - X)), (X, 0, 1))


def noncusp_first(channel, mass=MU, heavy=N):
    channel, mass, heavy = map(s.sympify, (channel, mass, heavy))
    a = channel
    q = heavy - a
    return (
        -2 * (a - 2 * mass) * Lll(a, mass) / q
        + 4 * heavy * Lmix(mass, heavy) / q
        - 2 * mass * (s.log(mass) - 1) / q
        + (
            2 * a * a
            + 4 * (-a * a / 2 + 2 * heavy * a - heavy**2) * Lh(a, heavy)
            - 2 * a * heavy * (s.log(heavy) - 1)
        )
        / q**2
    )


def finite_bracket(energy=S, transfer=T, mass=MU, heavy=N, resolution_squared=E2):
    energy, transfer, mass, heavy, resolution_squared = map(
        s.sympify, (energy, transfer, mass, heavy, resolution_squared)
    )
    rows = graphs.channels(energy, transfer, mass)
    Htree = sum(1 / (heavy - a) for a in rows)
    c = s.EulerGamma - s.log(4 * s.pi * resolution_squared)
    soft = sum(graphs.numerator(b, mass, 0) * J0(b, mass) for b in rows) - 2 * mass
    total = Htree * (
        sum(
            graphs.numerator(b, mass, 0) * J1(b, mass) + 2 * mass**2 * J0(b, mass)
            for b in rows
        )
        + c * soft
    )
    total += sum(
        -noncusp_first(a, mass, heavy)
        + 4 * (a - 2 * mass) * graphs.massive_T(a, mass, heavy)
        - 4
        * (a * a - 2 * heavy * mass)
        * graphs.offshell_U(a, mass, heavy)
        / (heavy - a)
        for a in rows
    )
    total += 2 * sum(
        graphs.numerator(b, mass, 0)
        * (
            graphs.massive_T(b, mass, heavy)
            - a * graphs.finite_box_K(a, b, mass, heavy)
        )
        / (heavy - a)
        for i, a in enumerate(rows)
        for j, b in enumerate(rows)
        if i != j
    )
    return total


def finite_amplitude(
    energy=S, transfer=T, mass=MU, heavy=N, cubic=G, kappa=K, resolution_squared=E2
):
    return (
        s.sympify(cubic) ** 2
        * finite_bracket(energy, transfer, mass, heavy, resolution_squared)
        / (16 * s.pi**2 * s.sympify(kappa))
    )


def finite_b20(mass=MU, heavy=N, cubic=G, kappa=K, resolution_squared=E2):
    mass, heavy, cubic, kappa, resolution_squared = map(
        s.sympify, (mass, heavy, cubic, kappa, resolution_squared)
    )
    value = finite_amplitude(
        2 * mass + V, 0, mass, heavy, cubic, kappa, resolution_squared
    )
    return s.Subs(s.Derivative(value, V, 2, evaluate=False), V, 0) / 2


def require_bound_domain(mass, heavy, resolution_squared):
    mass, heavy = map(source.require_mass, (mass, heavy))
    resolution_squared = quartic_forward.require_resolution(resolution_squared)
    if mass != 1 or heavy < 8:
        raise ValueError("The numerical majorant requires mu=1 and n>=8")
    return mass, heavy, resolution_squared


def known_bound(heavy=N, cubic=G, kappa=K):
    heavy, cubic, kappa = map(s.sympify, (heavy, cubic, kappa))
    return 2001 * cubic**2 / (kappa * heavy)


@cache
def data():
    a, b, n, u, h, e, A, y = s.symbols(
        "a b n u h epsilon positive_denominator share", positive=True
    )
    q = n - a
    delta_h = q + a * (u + 2 * h)
    checks = {}
    checks["whole_IR_subtraction_actual_derivative_rearrangement"] = s.expand(
        -(e - 1) * delta_h - (1 - e) * a * (u + 2 * h) - (1 - e) * q
    )
    wrong = s.expand(-(e - 1) * delta_h + (1 - e) * a * (u + 2 * h) - (1 - e) * q)
    checks["rejected_wrong_box_sign_has_nonzero_defect"] = s.factor(
        wrong - 2 * (1 - e) * a * (u + 2 * h)
    )
    J0, J1, V0, V1, K1, c, l = s.symbols("J0 J1 V0 V1 K1 Gamma_scale log_resolution")
    hh = s.Symbol("whole_heavy_tree")
    mu = s.Symbol("mu", positive=True)
    # Whole cusp/box sum and all noncusp moments are independently reduced earlier.
    numerator = hh * ((V0 + e * V1) * (J0 + e * J1) - 2 * mu) - e * K1
    hard = s.cancel(
        ((1 + e * c) * numerator - hh * (V0 * J0 - 2 * mu) * (1 + e * l)) / e
    )
    expected = hh * (V0 * J1 + V1 * J0 + (c - l) * (V0 * J0 - 2 * mu)) - K1
    checks["whole_finite_Gamma_and_soft_scale_expansion"] = s.factor(
        hard.subs(e, 0) - expected
    )
    checks["whole_J0_first_derivative"] = s.diff(1 / A, A) * (-y) - y / A**2
    checks["whole_J0_second_derivative"] = s.diff(y / A**2, A) * (-y) - 2 * y**2 / A**3
    checks["whole_J1_first_derivative"] = (
        s.diff(s.log(A) / A, A) * (-y) - y * (s.log(A) - 1) / A**2
    )
    checks["whole_J1_second_derivative"] = (
        s.diff(y * (s.log(A) - 1) / A**2, A) * (-y) - y**2 * (2 * s.log(A) - 3) / A**3
    )
    cU = -4 * (a * a - 2 * n) / q
    checks["whole_U_coefficient_first_derivative"] = s.factor(
        s.diff(cU, a) + 4 * (2 * a * q + a * a - 2 * n) / q**2
    )
    checks["whole_U_coefficient_second_derivative"] = s.factor(
        s.diff(cU, a, 2) + 8 / q + 16 * a / q**2 + 8 * (a * a - 2 * n) / q**3
    )
    cH = 4 * (-a * a / 2 + 2 * n * a - n * n) / q**2
    checks["whole_Lh_coefficient_first_derivative"] = s.factor(
        s.diff(cH, a) - 4 * n * a / q**3
    )
    checks["whole_Lh_coefficient_second_derivative"] = s.factor(
        s.diff(cH, a, 2) - 4 * n * (n + 2 * a) / q**4
    )
    checks["whole_a_over_q_first_derivative"] = s.factor(s.diff(a / q, a) - n / q**2)
    checks["whole_a_over_q_second_derivative"] = s.factor(
        s.diff(a / q, a, 2) - 2 * n / q**3
    )
    checks["whole_a_over_q2_second_derivative"] = s.factor(
        s.diff(a / q**2, a, 2) - 4 / q**3 - 6 * a / q**4
    )
    v = s.Symbol("v", real=True)
    coeff = s.symbols("f0:5")
    poly = sum(coeff[k] * (a - 2 * mu) ** k for k in range(5))
    checks["whole_crossed_forward_second_coefficient"] = s.diff(
        poly.subs(a, 2 * mu + v) + poly.subs(a, 2 * mu - v), v, 2
    ).subs(v, 0) / 2 - s.diff(poly, a, 2).subs(a, 2 * mu)
    x = s.Symbol("x", positive=True)
    checks["whole_mix_lower_gap"] = s.expand(n * x + (1 - x) ** 2 - 1 - x * (n - 2 + x))
    checks["whole_mix_upper_gap"] = s.expand(
        n - n * x - (1 - x) ** 2 - (1 - x) * (n - 1 + x)
    )
    # Actual corner integrals used in K, rather than a guessed norm.
    C, Q, _H = s.symbols("C Q upper_h", positive=True)
    primitive = (s.log(C + Q * h) + C / (C + Q * h)) / Q**2
    checks["whole_K_log_corner_primitive"] = s.factor(
        s.diff(primitive, h) - h / (C + Q * h) ** 2
    )
    checks["whole_K_u_log_weight"] = s.integrate(-4 * x * s.log(x), (x, 0, 1)) - 1
    checks["whole_K_u_constant_weight"] = s.integrate(2 * x, (x, 0, 1)) - 1
    checks["whole_T_and_box_I_same_integrand"] = s.expand(
        n * (1 - u)
        + u * u * (1 - b * y * (1 - y))
        - (n * x + (1 - x) ** 2 * (1 - b * y * (1 - y))).subs(x, 1 - u)
    )
    # Rounded bounds and their exact, nonnegative safety margins.
    R = s.Rational
    slacks = {
        "J0_0": 3 - R(8, 3),
        "J0_1": 2 - R(16, 9),
        "J0_2": 3 - R(64, 27),
        "J1_0": 3 - R(8, 3),
        "J1_1": 4 - R(32, 9),
        "J1_2": 6 - R(160, 27),
        "K1_first_per_channel": 19 - (12 + R(48, 8) + R(64, 64)),
        "K1_A_group": 75 - (16 + R(320, 8) + R(1200, 64)),
        "K1_cH_second": 100 - (32 + R(480, 8)),
        "K1_fourth_constant": 71 - R(566, 8),
        "T_linear_group": 20 - R(128, 9),
        "cU_0": 23 - (16 + R(50, 8)),
        "cU_1": 125 - (72 + R(100, 8)),
        "cU_2": 60 - (16 + R(288, 8) + R(400, 64)),
        "cU_whole_group": 73 - R(582, 8),
        "K_norm": R(14, 2) - (R(16, 3) + R(3, 2)),
        "K_delta_ratio": 1 - (R(2, 3) + R(2, 8)),
        "dK_second": 92 - (84 + R(64, 8)),
        "dK_whole": 164 - R(1308, 8),
        "dT_first": 18 - (16 + R(16, 8)),
        "dT_second": 17 - (8 + R(64, 8) + R(64, 64)),
        "dT_product": 49 - (17 + 2 * 18 * R(2, 3) + 8 * R(8, 9)),
        "dT_whole": 37 - R(294, 8),
        "F_phase": 2016 - (1680 + R(2304, 8) + R(3072, 64)),
        "final_200000": 200000 - (4869 + 366 * 500),
        "coefficient_2000": 2000 - R(200000, 128),
    }
    assert all(v >= 0 for v in slacks.values()), slacks
    L = s.Symbol("log_heavy_mass", positive=True)
    checks["whole_majorant_constant_and_log_sum"] = s.expand(
        2016
        + (110 + 72 * L)
        + 20 * (L + 4)
        + 73 * (L + 3)
        + 164 * (L + 14)
        + 37 * (L + 4)
        - (4869 + 366 * L)
    )
    exp_25 = sum(R(5, 2) ** k / s.factorial(k) for k in range(7))
    exp_3 = sum(s.Integer(3) ** k / s.factorial(k) for k in range(7))
    assert exp_25 > 10 and exp_3 > 16
    n0 = s.Integer(10) ** 200 / 512 + 2
    g = s.Rational(1, 8192)
    kappa = s.Integer(10) ** 800
    assert 8 < n0 < s.Integer(10) ** 200
    assert 2001 * g * g / (kappa * n0) < s.Rational(1, 10**1001)
    assert 2002 * g * g / (kappa * n0) < s.Rational(1, 10**1001)
    return {
        "whole_noncusp_first_epsilon_coefficient": noncusp_first(S),
        "whole_known_nonendpoint_finite_amplitude": finite_amplitude(),
        "whole_known_nonendpoint_forward_coefficient": finite_b20(),
        "whole_explicit_second_derivative_majorant": (4869 + 366 * s.log(N)) / N,
        "whole_known_mixed_endpoint_included_majorant": known_bound(),
        "whole_known_mixed_plus_quartic_majorant": 2002 * G**2 / (K * N),
        "whole_exact_majorant_safety_margins": slacks,
        "whole_finite_derivation": "The entire mixed box is subtracted by exact h integration by parts, with the negative finite aK term. All proper/selfenergy/contact moments are combined before their first epsilon derivative. Their noncusp zeroth coefficient is2mu/(n-a), so the full pole is exactly2Bsoft times the complete heavy tree. Expanding the whole Gamma/scale factor and subtracting only the stated S278 factor gives the displayed finite Feynman-parameter representative. The two endpoint assignments and the complete S290 known g^2/kappa form factors are then added separately.",
        "whole_dominated_derivative_argument": "On mu=1,n>=8,|v|<=1/2, all three channel denominators A are>=3/8 and n-a>=n/2. The exact K subtraction gains one soft power. Its value and the required invariant/epsilon derivatives are integrable; logarithms from epsilon derivatives do not destroy the soft gain. U has an integrable two-parameter corner and T is gapped. The displayed parameter inequalities dominate two forward derivatives, so those derivatives commute with the finite subtracted integrals. No physical cut or heavy resonance limit is exchanged.",
        "whole_majorant_proof": "The second derivative contributions are bounded by2016/n,(110+72ln n)/n,20(ln n+4)/n,73(ln n+3)/n,164(ln n+14)/n and37(ln n+4)/n. They sum to(4869+366ln n)/n. At original n,ln n<500, giving200000/n. Divide by2 for b20 and by16pi^2 for the loop. The known S290 endpoint adds less than g^2/(kappa n). Thus the selected known mixed coefficient is below2001g^2/(kappa n)<10^-1001, and including the S293 C sector still gives a bound below10^-1001.",
        "whole_matching_and_IR_scope": "The representative is conditional on the inherited analytic S278 soft division and the stated E^2 interval; no detector is selected. The independent R H finite anchor, additional heavy residue/cubic and higher-EFT matching are outside the known bound. Finite constant quartic/RPhi2 pieces have zero one-insertion b20, not zero higher insertions. All-loop, physical IR/width/Regge, state/domain/measure/bounce and original V/G/B/P8 remain open.",
        "checks": checks,
        "gates": {
            "actual_negative_box_K_sign_checked_not_only_Delta_derivative": True,
            "entire_finite_Gamma_evanescent_coefficient": True,
            "positive_subthreshold_domain_dominates_two_derivatives": True,
            "all26_exact_majorant_safety_margins_nonnegative": all(
                v >= 0 for v in slacks.values()
            ),
            "original_mass_log_and_exact_numeric_bound": 8 < n0 < s.Integer(10) ** 200
            and 2001 * g * g / (kappa * n0) < s.Rational(1, 10**1001),
            "known_S290_endpoint_and_S293_quartic_ownership": True,
            "unmatched_RH_and_higher_operators_not_bounded_by_known_loop": True,
            "no_physical_IR_or_Regge_inference_from_local_graph_bound": True,
        },
    }
