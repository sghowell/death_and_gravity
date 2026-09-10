"""Mixed heavy/spectral bubble with the complete outer on-shell forest."""

from functools import cache

import sympy as sp
from p8_vacuum_light_pole import kernel as original


@cache
def data():
    x, s, u, M, g, Q = sp.symbols("x s spectral_mass_squared heavy_mass_squared g Q")
    A = x * (1 - x)
    Delta = x * u + (1 - x) * M - A * s
    Delta1 = Delta.subs(s, 1)
    b = A / Delta1
    K = -sp.log(Delta)
    raw = K - K.subs(s, 1) - (s - 1) * sp.diff(K, s).subs(s, 1)
    R = -sp.log(1 - (s - 1) * b) - (s - 1) * b
    d = sp.symbols("positive_spacelike_displacement", positive=True)
    real_R = d * b - sp.log(1 + d * b)
    old = original.data()
    frozen_delta = old["Delta"].subs(
        {
            old["Feynman_parameter"]: 1 - x,
            old["Minkowski_invariant"]: s,
            old["heavy_mass_squared"]: M,
        },
        simultaneous=True,
    )
    checks = {
        "same_frozen_mass_one_light_bubble_with_parameter_reflection": sp.factor(
            frozen_delta - Delta.subs(u, 1)
        ),
        "same_frozen_mixed_bubble_prefactor": sp.factor(
            old["self_energy_prefactor"] * (16 * sp.pi**2) - old["cubic_squared"]
        ),
        "whole_OS_remainder_derivative_agreement": sp.factor(sp.diff(raw - R, s)),
        "whole_OS_remainder_anchor_agreement": sp.simplify((raw - R).subs(s, 1)),
        "outer_pole_mass_anchor": R.subs(s, 1),
        "outer_unit_residue_anchor": sp.diff(R, s).subs(s, 1),
        "outer_first_derivative_kernel": sp.factor(sp.diff(K, s) - A / Delta),
        "outer_second_derivative_kernel": sp.factor(
            sp.diff(K, s, 2) - A * A / Delta**2
        ),
        "halfplane_gap_margin": sp.expand(
            Delta - x * (u - 4) - (1 - x) * M - A * (4 - s) - 4 * x * x
        ),
        "spacelike_parameter_gap": sp.factor(
            Delta.subs(s, 1 - d) - Delta1 * (1 + d * b)
        ),
        "spacelike_remainder_from_same_OS_kernel": sp.simplify(
            R.subs(s, 1 - d) - real_R
        ),
        "spacelike_positive_second_derivative": sp.factor(
            sp.diff(real_R, d, 2) - b * b / (1 + d * b) ** 2
        ),
        "spacelike_remainder_below_slope_times_displacement": sp.factor(
            d * b - real_R - sp.log(1 + d * b)
        ),
    }
    return {
        "symbols": {"x": x, "s": s, "u": u, "M": M, "g": g, "Q": Q},
        "mixed_mass_parameter_denominator": Delta,
        "on_shell_parameter_gap": Delta1,
        "whole_outer_OS_parameter_kernel": raw,
        "anchored_analytic_OS_kernel": R,
        "positive_spacelike_OS_kernel": real_R,
        "finite_outer_slope_integrand": b,
        "spectral_two_point_family": "Pi_spec,R(s)=(g/Q) integral_T^infinity w(u) integral_0^1 R(s,x,u) dx du. w is the whole first fermion insertion measure of S6.135.",
        "inverse_sign": "The continued Euclidean inverse receives -Pi_spec,R, as for the original mixed H-Phi bubble. There is one spectral light line, not two.",
        "domain": "M>=1, mF>=2, u>=T=4mF^2; Re(s)<=4 has a uniform positive parameter gap. H is a perturbative heavy mass reference, not a declared exact stable particle.",
        "scope": "On-shell-renormalized first covariance-insertion quadratic family. Other two-loop primitive graphs and reference conversions are not included.",
        "checks": checks,
    }
