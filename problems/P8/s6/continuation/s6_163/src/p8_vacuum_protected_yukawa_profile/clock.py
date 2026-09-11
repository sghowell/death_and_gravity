"""Conditional mass-derivative bounds, not a transported bounce state."""

from functools import cache

import sympy as s
from p8_offshell_vacuum import jets as original_jets

from . import profile


def ratios(mass, yukawa_upper, radius, velocity, acceleration, eta=profile.ETA):
    m, y, r, v, a, e = map(
        profile.exact, (mass, yukawa_upper, radius, velocity, acceleration, eta)
    )
    if min(m, y, r) <= 0 or min(v, a) < 0 or not 0 < e < 1:
        raise ValueError(
            "Require positive mass/coupling/radius and nonnegative argument derivative bounds"
        )
    if not y * r < e * m:
        raise ValueError(
            "The supplied first-power Yukawa cap does not prove this mass floor"
        )
    floor = (1 - e) * m
    return {
        "first_mass_derivative_over_mass_squared_upper": y * v / floor**2,
        "second_mass_derivative_over_mass_cubed_upper": y
        * (a + 9 * v * v / r)
        / floor**3,
    }


@cache
def data():
    R, v, a, y = s.symbols("R v a y", positive=True)
    z = s.Function("z")
    t = s.Symbol("t", real=True)
    f = profile.expression(z(t), R)
    chain = s.diff(f, t, 2)
    second = profile.data()["argument_second_derivative"].subs(
        {profile.Z: z(t), profile.R: R}
    )
    first = profile.data()["argument_first_derivative"].subs(
        {profile.Z: z(t), profile.R: R}
    )
    phi, kappa = s.symbols("Phi kappa", positive=True)
    lam, gamma, c = s.symbols("quartic_lambda quartic_gamma redundant_c", real=True)
    old = original_jets.data()
    substitutions = {
        symbol: (
            phi if sum(index) == 0 else s.sqrt(kappa) if index == (1, 0, 0, 0) else 0
        )
        for symbol, index in original_jets.BY_SYMBOL.items()
    }
    flat_map = s.expand(
        phi + old["cubic_field_redefinition"].subs(substitutions, simultaneous=True)
    )
    linear = 1 + 2 * lam * kappa - 2 * gamma * kappa
    cubic = lam - s.Rational(2, 3) * gamma - c / 6
    return {
        "conditional_chain_bounds": {"first": y * v, "second": y * (a + 9 * v * v / R)},
        "flat_linear_clock_extrapolation_of_literal_cubic_map": flat_map,
        "flat_linear_clock_extrapolated_map_linear_coefficient": linear,
        "diagnostic_curve": "The numerical specialization z(t)=sqrt(kappa)t is a test curve for the saturation argument only. It is not asserted to equal the transported polynomial field F(Psi) on the actual bounce.",
        "scope": "These are pointwise first/second mass-derivative ratios conditional on the stated argument derivatives. They do not bound particle creation, a Hadamard-state remainder, metric time derivatives or the full Dirac propagator on a rolling geometry.",
        "checks": {
            "second_argument_chain_rule": s.simplify(
                chain - second * s.diff(z(t), t) ** 2 - first * s.diff(z(t), t, 2)
            ),
            "second_profile_numerator_coefficient": s.diff(
                profile.data()["argument_second_derivative"]
                * profile.R
                * (1 + (profile.Z / profile.R) ** 8) ** s.Rational(17, 8),
                profile.Z,
                7,
            ).subs(profile.Z, 0)
            + 9 * s.factorial(7) / profile.R**7,
            "real_second_derivative_cap_small_argument": 1**7 - 1,
            "real_second_derivative_cap_large_argument_exponent": 17 - 7 - 10,
            "literal_flat_linear_clock_map": s.expand(
                flat_map - linear * phi - cubic * phi**3
            ),
            "literal_flat_linear_clock_map_slope_at_zero": s.diff(flat_map, phi).subs(
                phi, 0
            )
            - linear,
        },
    }
