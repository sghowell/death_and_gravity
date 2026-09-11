"""Private complete curvature Euler tensors and quantitative jet majorants."""

from functools import cache
from itertools import product
from math import comb

import sympy as s

from . import metrics


@cache
def constants():
    G = metrics.G
    connection = [s.Rational(3, 2) * 2**j * G**2 for j in range(4)]
    riemann = [
        2 * connection[j + 1]
        + 2 * sum(comb(j, l) * connection[l] * connection[j - l] for l in range(j + 1))
        for j in range(3)
    ]
    Riem = 100 * G**4
    scalar = 400 * G**5
    hess_scalar = scalar + 2 * G**2 * scalar
    hess_ricci = (
        Riem
        + 2 * (3 * G**2) * Riem
        + 2 * (2 * G**2) * Riem
        + 3 * (2 * G**2) * (Riem + 2 * (2 * G**2) * Riem)
    )
    HessR = 1200 * G**7
    BoxR = 1200 * G**8
    BoxRic = 5000 * G**9
    Scalar = 400 * G**5
    Ric = 100 * G**4
    RicUp = 100 * G**6
    RiemLow = 100 * G**5
    Ric2 = 10000 * G**10
    eR = Ric + G * Scalar / 2
    eR2 = 2 * HessR + 2 * Scalar * Ric + G * (2 * BoxR + Scalar**2 / 2)
    eRic2 = HessR + BoxRic + 2 * RicUp * RiemLow + G * (BoxR / 2 + Ric2 / 2)
    finite = (
        s.Rational(5, 3) * 1000**2 * (300 * G**6)
        + (s.Rational(1, 30) + s.Rational(1, 15)) * 10**6 * G**11
    )
    final = 8 * 10**10 * G**14
    return {
        "Christoffel_time_jet_majorants": connection,
        "complete_Riemann_time_jet_majorants": riemann,
        "Riemann_uniform_display": Riem,
        "scalar_time_jet_display": scalar,
        "complete_scalar_covariant_second_before_rounding": hess_scalar,
        "complete_Ricci_covariant_second_before_rounding": hess_ricci,
        "Einstein_before_rounding": eR,
        "R_squared_Euler_before_rounding": eR2,
        "Ricci_squared_Euler_before_rounding": eRic2,
        "fixed_finite_covariant_Euler_before_rounding": finite,
        "complete_local_current_raw_parameter_display_before_rounding": final,
    }


def geometry(raw):
    """All homogeneous coordinate curvature jets and full covariant Euler tensors."""
    if len(raw) != 5 or any(g.shape != (4, 4) for g in raw):
        raise ValueError("Require five complete four-dimensional metric time jets")
    inverse = [raw[0].inv()]
    for n in range(1, 5):
        inverse.append(
            -inverse[0]
            * sum(
                (comb(n, j) * raw[j] * inverse[n - j] for j in range(1, n + 1)),
                s.zeros(4),
            )
        )
    gamma = []
    for n in range(4):
        out = {}
        for a, b, c in product(range(4), repeat=3):
            value = 0
            for j in range(4):
                for l in range(n + 1):
                    derivative = (
                        (raw[n - l + 1][j, c] if b == 0 else 0)
                        + (raw[n - l + 1][j, b] if c == 0 else 0)
                        - (raw[n - l + 1][b, c] if j == 0 else 0)
                    )
                    value += (
                        s.Rational(1, 2) * comb(n, l) * inverse[l][a, j] * derivative
                    )
            out[a, b, c] = s.cancel(value)
        gamma.append(out)
    Riem = []
    for n in range(3):
        out = {}
        for a, b, c, d in product(range(4), repeat=4):
            value = (gamma[n + 1][a, d, b] if c == 0 else 0) - (
                gamma[n + 1][a, c, b] if d == 0 else 0
            )
            value += sum(
                comb(n, l)
                * (
                    gamma[l][a, c, e] * gamma[n - l][e, d, b]
                    - gamma[l][a, d, e] * gamma[n - l][e, c, b]
                )
                for e in range(4)
                for l in range(n + 1)
            )
            out[a, b, c, d] = s.cancel(value)
        Riem.append(out)
    Ric = [
        s.Matrix(4, 4, lambda b, d, n=n: sum(Riem[n][a, b, a, d] for a in range(4)))
        for n in range(3)
    ]
    scalar = [
        s.cancel(
            sum(
                comb(n, l) * inverse[l][a, b] * Ric[n - l][a, b]
                for a, b in product(range(4), repeat=2)
                for l in range(n + 1)
            )
        )
        for n in range(3)
    ]
    HessR = s.Matrix(
        4,
        4,
        lambda b, a: (scalar[2] if a == b == 0 else 0) - gamma[0][0, b, a] * scalar[1],
    )
    BoxR = s.cancel(
        sum(inverse[0][a, b] * HessR[a, b] for a, b in product(range(4), repeat=2))
    )
    cov1 = {}
    cov1dot = {}
    for a, i, j in product(range(4), repeat=3):
        cov1[a, i, j] = (Ric[1][i, j] if a == 0 else 0) - sum(
            gamma[0][k, a, i] * Ric[0][k, j] + gamma[0][k, a, j] * Ric[0][i, k]
            for k in range(4)
        )
        cov1dot[a, i, j] = (Ric[2][i, j] if a == 0 else 0) - sum(
            gamma[1][k, a, i] * Ric[0][k, j]
            + gamma[1][k, a, j] * Ric[0][i, k]
            + gamma[0][k, a, i] * Ric[1][k, j]
            + gamma[0][k, a, j] * Ric[1][i, k]
            for k in range(4)
        )
    BoxRic = s.zeros(4)
    for b, a, i, j in product(range(4), repeat=4):
        second = (cov1dot[a, i, j] if b == 0 else 0) - sum(
            gamma[0][k, b, a] * cov1[k, i, j]
            + gamma[0][k, b, i] * cov1[a, k, j]
            + gamma[0][k, b, j] * cov1[a, i, k]
            for k in range(4)
        )
        BoxRic[i, j] += inverse[0][b, a] * second
    RicUp = inverse[0] * Ric[0] * inverse[0]
    mixed = s.Matrix(
        4,
        4,
        lambda i, j: sum(
            raw[0][i, a] * Riem[0][a, p, j, q] * RicUp[p, q]
            for a, p, q in product(range(4), repeat=3)
        ),
    )
    Ric2 = s.trace(RicUp * Ric[0])
    R = scalar[0]
    Einstein = Ric[0] - raw[0] * R / 2
    ER2 = 2 * HessR - 2 * R * Ric[0] + raw[0] * (-2 * BoxR + R * R / 2)
    ERic2 = HessR - BoxRic - 2 * mixed + raw[0] * (-BoxR / 2 + Ric2 / 2)
    return {
        "inverse": inverse,
        "Christoffel": gamma,
        "Riemann": Riem,
        "Ricci": Ric,
        "scalar": scalar,
        "scalar_Hessian": HessR,
        "scalar_box": BoxR,
        "Ricci_box": BoxRic.applyfunc(s.cancel),
        "Euler_R": (-Einstein).applyfunc(s.cancel),
        "Euler_R2": ER2.applyfunc(s.cancel),
        "Euler_Ricci2": ERic2.applyfunc(s.cancel),
    }


@cache
def data():
    c = constants()
    G = metrics.G
    r0, r1, r2, g0, g1 = s.symbols("T0 T1 T2 Gamma0 Gamma1", positive=True)
    rank = s.Integer(2)
    full = (
        r2 + rank * g1 * r0 + rank * g0 * r1 + (rank + 1) * g0 * (r1 + rank * g0 * r0)
    )
    checks = {
        "full_rank_two_covariant_derivative_count": s.expand(
            full - r2 - 2 * g1 * r0 - 5 * g0 * r1 - 6 * g0**2 * r0
        ),
        "canonical_fixed_scalar_curvature_density_reduction": s.Rational(-4, 1)
        * (s.Rational(1, 120))
        - s.Rational(-1, 30),
        "canonical_fixed_Ricci_density_reduction": s.Rational(-4, 1) * s.Rational(1, 60)
        - s.Rational(-1, 15),
        "raw_amplitude_second_factorial": s.factorial(2) - 2,
    }
    return {
        "tensor_norm": "All bounds use component-l1 factorial-weighted mixed time/amplitude jets. Tensor products and index contractions are bounded before contraction; homogeneous partial spatial derivatives do not remove covariant spatial connection terms.",
        "full_local_Euler": "Modulo compact four-dimensional Euler variation, the shear-dependent finite action is (5/3)m^2 R-R^2/30-Ricci^2/15. Use the complete covariant metric Euler tensors and then raise both metric indices for the detector pairing.",
        "cosmological_term": "The m^4 volume term has identically zero unimodular shear variation before any estimate, through all amplitude orders.",
        "constants": c,
        "checks": checks,
        "gates": {
            "full_Riemann_time_jets": all(
                v < c["Riemann_uniform_display"]
                for v in c["complete_Riemann_time_jet_majorants"]
            ),
            "scalar_covariant_second_display": c[
                "complete_scalar_covariant_second_before_rounding"
            ]
            < 1200 * G**7,
            "Ricci_covariant_second_display": c[
                "complete_Ricci_covariant_second_before_rounding"
            ]
            < 5000 * G**8,
            "Einstein_display": c["Einstein_before_rounding"] < 300 * G**6,
            "R_squared_Euler_display": c["R_squared_Euler_before_rounding"]
            < 10**6 * G**11,
            "Ricci_squared_Euler_display": c["Ricci_squared_Euler_before_rounding"]
            < 10**6 * G**11,
            "fixed_mass_complete_finite_Euler_display": c[
                "fixed_finite_covariant_Euler_before_rounding"
            ]
            < 10**10 * G**11,
            "all_local_raw_amplitude_derivatives_below_display": c[
                "complete_local_current_raw_parameter_display_before_rounding"
            ]
            < s.Integer(10) ** 86,
        },
    }
