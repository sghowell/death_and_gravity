"""Fraction-only full sampling costs and complete-control reconstruction."""

from fractions import Fraction as F
from math import comb, factorial


def norms(length):
    return F(13, 35) * length, F(6, 5) / length, F(12) / length**3


def interval(length, lo, hi, ds):
    n0, n1, n2 = norms(length)
    d1, d2, d3, d4 = ds
    P = 3 * (hi * hi * n2 + 4 * d1 * d1 * n1 + 4 * (d1 * d1 / lo + d2) ** 2 * n0)
    M = 2 * (hi * hi * n1 + 4 * d1 * d1 * n0) / lo**4
    return {
        "P_norm_majorant": P,
        "mixed_norm_majorant": M,
        "quantum_cost": F(5, 9) * P + M,
        "state_cost": 2 * n1 + 8 * d1 * d1 * n0 / lo**2 + 3 * n0 / lo**4,
        "reference_zero_cost": 4 * d1 * d1 * d2 / lo * n0 / 360,
        "reference_beta_cost": (
            48 * d1 * d1 * d2 / lo + 84 * d2 * d2 + 72 * d1 * d3 + 12 * hi * d4
        )
        * n0
        / 360,
        "sampler_mass": n0,
    }


def cost_replay():
    past = interval(F(1, 1000), F(9, 10), F(2), list(map(F, (3, 16, 1000, 2000000))))
    future = interval(
        F(2), F(1, 2), F(2), list(map(F, (50, 10000, 2000000, 400000000)))
    )
    quantum = past["quantum_cost"] + future["quantum_cost"]
    c0 = quantum + past["reference_zero_cost"] + future["reference_zero_cost"]
    cbeta = past["reference_beta_cost"] + future["reference_beta_cost"]
    state = past["state_cost"] + future["state_cost"]
    mass = past["sampler_mass"] + future["sampler_mass"]
    ricci = F(320, 9) * past["sampler_mass"]
    kinetic = F(6, 5)
    upper = kinetic + ricci + max(c0, cbeta) / 10**12 + state / 10**6 + mass / 10
    return {
        "past": past,
        "future": future,
        "quantum_cost": quantum,
        "C0": c0,
        "Cbeta": cbeta,
        "Cstate": state,
        "Csource": mass,
        "past_Ricci_upper_cost": ricci,
        "future_index_kinetic_cost": kinetic,
        "worst_gate_index_upper": upper,
        "strict_index_margin": F(2) - upper,
    }


def control_replay():
    k = F(100)
    r = F(1, 1000)
    q = [
        F(0),
        F(2),
        2 * k + 4,
        k * k + 4 * k + F(40, 3),
        k**3 / 3 + 2 * k * k + F(40, 3) * k + F(160, 3),
    ]
    S = sum((q[n] * F(factorial(n)) / k ** (n - 1) for n in range(1, 5)), F(0))
    future = {j: 2**j * k ** (j - 1) * S for j in range(1, 5)}
    past = {}
    # (d/dx-k)^j Q gives the derivative polynomial multiplying exp(-kx).
    for j in range(1, 5):
        coeff = [
            sum(
                (
                    F(comb(j, m))
                    * (-k) ** (j - m)
                    * F(factorial(n + m), factorial(n))
                    * q[n + m]
                    for m in range(j + 1)
                    if n + m < len(q)
                ),
                F(0),
            )
            for n in range(len(q))
        ]
        past[j] = F(10, 9) * sum((abs(c) * r**n for n, c in enumerate(coeff)), F(0))
    return {
        "future_S": S,
        "future_A_lower": 1 - S / k,
        "future_derivative_bounds": future,
        "past_absolute_A_minus_one_bound": F(10, 9)
        * sum((abs(c) * r**n for n, c in enumerate(q)), F(0)),
        "past_derivative_bounds": past,
        "past_derivative_margins": {
            j: F((3, 16, 1000, 2000000)[j - 1]) - past[j] for j in range(1, 5)
        },
        "future_derivative_margins": {
            j: F((50, 10000, 2000000, 400000000)[j - 1]) - future[j]
            for j in range(1, 5)
        },
    }


def replay():
    return {"costs": cost_replay(), "complete_geometry_control": control_replay()}
