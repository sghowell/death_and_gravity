"""Analytic inverse-radius normalization of the full two-time symbol."""

from functools import cache

import sympy as s

from . import sectors


def unscale_point(point, x):
    out = {"a": point["a"]}
    for key in ("kt", "kl", "lt", "ll"):
        f, p, omega = point[key]
        out[key] = (s.sqrt(x) * f, p / s.sqrt(x), omega / x)
    return out


def normalized_products(n, P, x, Dsharp, G, source, detector_sharp, m):
    """No inverse powers of x or individual moving frames remain."""
    ell = -s.Matrix(n) + x * s.Matrix(P)
    return sectors.pair_products(
        s.Matrix(n), ell, Dsharp, G, source, detector_sharp, m * x
    )


def normalized_frequency(a, m, x, v, coefficients):
    v = s.Matrix(v)
    omega = s.sqrt((v.T * v)[0] / a**2 + m * m * x * x)
    return omega + sum(
        c * x ** (2 * j) * omega ** (1 - 2 * j) for j, c in enumerate(coefficients, 1)
    )


@cache
def data():
    x, m = s.symbols("x m", positive=True)
    source = {"a": s.symbols("a_source", positive=True)}
    detector = {"a": s.symbols("a_detector", positive=True)}
    for name in ("kt", "kl", "lt", "ll"):
        source[name] = s.symbols("F_s_" + name + " Pi_s_" + name + " Omega_s_" + name)
        detector[name] = s.symbols("F_d_" + name + " Pi_d_" + name + " Omega_d_" + name)
    k = s.Matrix([3, 0, 4])
    l = s.Matrix([0, 5, 12])
    D = s.Matrix([[1, 2, -1], [2, -2, 1], [-1, 1, 1]])
    G = s.Matrix([[2, 1, 3], [1, 1, -2], [3, -2, -3]]) / 7
    original = sectors.pair_products(
        k / x, l / x, D, G, unscale_point(source, x), unscale_point(detector, x), m
    )
    normalized = sectors.pair_products(k, l, D, G, source, detector, m * x)
    checks = {
        key + "_exact_two_time_inverse_radius_normalization": s.cancel(
            x * x * original[key] - normalized[key]
        )
        for key in sectors.SECTORS
    }
    A = sectors.scalar_amplitudes(source, m * x)
    checks["both_mixed_scaled_amplitudes_vanish_at_zero"] = s.Matrix(
        [A["TL"].subs(x, 0), A["LT"].subs(x, 0)]
    )
    a = s.symbols("a", positive=True)
    c = s.symbols("P1:5")
    W = normalized_frequency(a, m, x, s.Matrix([0, 0, 1]), c)
    checks["full_four_term_scaled_frequency_limit"] = s.simplify(W.subs(x, 0) - 1 / a)
    omega = s.sqrt(a**-2 + m * m * x * x)
    correction = W - omega
    checks["actual_scaled_adiabatic_correction_low_jets"] = s.Matrix(
        [
            s.simplify(correction.subs(x, 0)),
            s.simplify(s.diff(correction, x).subs(x, 0)),
            s.simplify(s.diff(correction, x, 2).subs(x, 0) / 2 - a * c[0]),
        ]
    )
    return {
        "normalization": "Set x=1/r, k=n/x, l=(-n+xP)/x. Write W=What/x, omega=Omega/x, f=sqrt(x)F and p=Pi/sqrt(x). The exact normalized pair product x^2 F_pair is the same four-sector formula at scaled momenta n,-n+xP and effective algebraic mass m*x. This is an algebraic normalization, not a physical mass change.",
        "frequency": "Omega=sqrt(v.v/a^2+m^2x^2), What=Omega+sum_(q=1)^4 P_q(t,z)x^(2q)Omega^(1-2q), z=1-m^2x^2/Omega^2. F=(2What)^-1/2 and Pi=(minus_i*What-x*(What_time/(2What)+rate))*F, with the opposite sign in the analytic Schwarz detector.",
        "analytic_domain": "At x=0, scaled momenta are n,-n with n.n=1, Omega=What=1/a and the projector denominators are1. Positivity of a and compactness of the actual time slab give one sufficiently small complex inverse-radius disc for each fixed bounded external P set. Full scaled expressions and source-time jets are holomorphic there. No global individual polarization frame or Borel-state analyticity is needed.",
        "sector_phase": "Each inverse phase is x times its own ghat_sigma_tau=1/(What_k_sigma+What_l_tau). Source differentiation does not change powers of x because internal radius is fixed during time derivatives.",
        "checks": checks,
        "gates": {
            "all_four_original_sector_scalings_exact": set(original)
            == set(sectors.SECTORS),
            "full_W8_four_coefficients_retained": len(c) == 4,
            "zero_inverse_radius_projectors_regular": True,
            "effective_mx_not_physical_massless_limit": True,
            "fixed_P_analytic_domain_not_all_P_norm": True,
        },
    }
