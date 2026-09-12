"""Three-subtracted nonlocal cut remainder, with local matching left explicit."""

from functools import cache

import sympy as s

from . import projectors as pro


def moment(spin, order):
    pro.above_threshold_density(spin)
    if isinstance(order, bool) or not isinstance(order, (int, s.Integer)) or order < 3:
        raise ValueError("An exact moment order of at least3 is required")
    return _moment(int(spin), int(order))


@cache
def _moment(spin, order):
    v = s.Symbol("v", nonnegative=True)
    m = pro.MASS
    sigma = 4 * m * m / (1 - v * v)
    polynomial = (
        (13 * sigma * sigma + 56 * m * m * sigma + 48 * m**4) / 3840
        if spin == 2
        else (sigma * sigma - 4 * m * m * sigma + 12 * m**4) / 384
    )
    transformed = s.cancel(
        v * polynomial * (8 * m * m * v / (1 - v * v) ** 2) / sigma ** (order + 1)
    )
    return s.integrate(transformed, (v, 0, 1)) / s.pi**2


@cache
def data():
    sigma, z = s.symbols("sigma z")
    m = pro.MASS
    checks = {
        "complete_three_subtraction_identity": s.factor(
            1 / (sigma - z)
            - 1 / sigma
            - z / sigma**2
            - z * z / sigma**3
            - z**3 / (sigma**3 * (sigma - z))
        ),
        "spin2_first_convergent_moment": s.simplify(
            moment(2, 3) - 3 / (3584 * s.pi**2 * m * m)
        ),
        "spin0_first_convergent_moment": s.simplify(
            moment(0, 3) - 3 / (8960 * s.pi**2 * m * m)
        ),
        "positive_moment_ratio": s.simplify(
            moment(2, 3) / moment(0, 3) - s.Rational(5, 2)
        ),
        "spectral_tail_inverse_lower_limit": s.integrate(
            sigma**-2, (sigma, s.Symbol("Lambda", positive=True), s.oo)
        )
        - 1 / s.Symbol("Lambda", positive=True),
    }
    return {
        "definition": "D_i(z)=z^3 integral_(4m^2)^infinity rho_i(sigma)/[sigma^3(sigma-z)]dsigma. This is the three-subtracted dispersive PART; the physical polynomial c0+c1 z+c2 z^2 and local tensor/contact matching are not chosen here.",
        "analyticity": "The positive density grows as sigma^2. Three subtractions give a locally dominated integral analytic on C minus[4m^2,infinity). Two subtractions do not give an absolutely integrable ultraviolet majorant.",
        "controlled_low_energy": "For |z|<=2m^2, |D_i(z)|<=2|z|^3 I_i, I_i=integral rho_i/sigma^4. Also |D_i-z^3 I_i|<=2|z|^4 J_i, J_i=integral rho_i/sigma^5. These are convergent remainder bounds, not a chosen physical subtraction polynomial.",
        "moments": {
            "spin2_I": moment(2, 3),
            "spin0_I": moment(0, 3),
            "spin2_J": moment(2, 4),
            "spin0_J": moment(0, 4),
        },
        "UV_tail": "For a computational spectral lower limit Lambda>=4m^2 and |z|<=2m^2, the discarded D2 tail is<=2|z|^3/(128pi^2 Lambda) and D0 tail<=2|z|^3/(384pi^2 Lambda). Lambda is a squared spectral mass and not a physical cutoff.",
        "boundary": "The absorptive cut cannot fix finite local counterterms, full Ward/contact matching, a curved-background retarded kernel, a feedback inverse or finite-gravity Regge information.",
        "checks": checks,
        "gates": {
            "both_first_moments_positive": moment(2, 3) > 0 and moment(0, 3) > 0,
            "both_next_moments_positive": moment(2, 4) > 0 and moment(0, 4) > 0,
            "threshold_stays_massive": m.is_positive is True,
        },
    }
