"""Complete analytic canonical source on a uniform finite-jet polydisc."""

from functools import cache

import sympy as s
from p8_vacuum_canonical_affine_decoupling import family
from p8_vacuum_canonical_affine_decoupling import source as old

N, K0, MASS = family.N, family.K0, s.Integer(1000)
K = s.Symbol("kappa", positive=True)
M = 20000 * N / (K * K0)
C = 30 * M


@cache
def data():
    n, k0 = s.symbols("n kappa0", positive=True)
    rho = 1 / (16 * n)
    V = 10240 * n / k0
    B = 307200 * n / k0**2
    Ct = 943718400 * n**2 / k0**3
    normalized = s.factor((V + B + Ct) / (n / k0))
    endpoint = s.factor(normalized.subs(k0, 512 * n))
    pm, Y, Q, Z, a, ay, ap, H, R = s.symbols(
        "Phi_mu Y BoxPhi ZPhi regular_a a_Y a_Phi Href Rk", real=True
    )
    kk = s.Symbol("kappa", positive=True)
    U = (
        pm * a * (Y * Q - Z)
        - 3 * pm * Y**2 * a * H / s.sqrt(kk)
        + 3 * pm * Y**2 * a * (Y**2 * ap / 4 + (2 * a + Y * ay) * Z / 2) / (kk * R)
    )
    prior = old.data()["regular_U_mu"]
    mapping = {str(x): x for x in prior.free_symbols}
    table = {
        "Phi_mu": pm,
        "Y": Y,
        "Box_Phi": Q,
        "Z_Phi": Z,
        "regular_a": a,
        "a_Y": ay,
        "a_Phi": ap,
        "Href": H,
        "kappa": kk,
    }
    actual = prior.subs(
        {mapping[key]: value for key, value in table.items()}, simultaneous=True
    )
    checks = {
        "literal_full_canonical_source": s.factor(
            actual - U.subs(R, 1 + Y * Y * a / kk)
        ),
        "outer_regular_A_product": 20 * n - 5 * n * 2 * 2,
        "inner_A_u_Cauchy": 160 * n - 20 * n / s.Rational(1, 8),
        "inner_A_X_Cauchy": 640 * n**2 - 20 * n / (rho / 2),
        "canonical_V_bound": V - 2 * (20 * n / k0) * (16 * 8 + 128),
        "canonical_B_over_sqrt_kappa_bound": B
        - 3 * 2 * 256 * (20 * n / k0) * (10 / s.sqrt(k0)) / s.sqrt(k0),
        "canonical_C_over_kappa_R_bound": Ct
        - 2 * (3 * 2 * 256 * 20 * n / k0) * (15360 * n / k0) / k0,
        "normalized_polydisc_source": normalized
        - 10240
        - 307200 / k0
        - 943718400 * n / k0**2,
        "minimum_hierarchy_source_enclosure": endpoint - 10240 - 4200 / n,
        "fifteen_independent_scalar_jets": 1 + 4 + 10 - 15,
        "four_component_first_derivative_sum": 30 - 2 * 15,
        "all_spatial_source_derivatives": 60 - 2 * 30,
        "all_spatial_Frechet_derivatives": 1860 - 2 * 2 * (15 + 2 * 15**2),
        "Frechet_spatial_to_undifferentiated_ratio": s.Rational(1860, 30) - 62,
    }
    return {
        "full_source": "S_mu=U_mu/kappa from literal S6.177 regular all-Y source, including Href and every R-dependent term",
        "same_action": "S6.182 changes only fixed F; R, complete affine source and canonical vector are unchanged",
        "canonical_class": "On a time interval of length at most one, real smooth Phi is prepared zero initially and spatially Schwartz, with every coordinate jet through order3 bounded pointwise by1. No time Fourier band or finite source harmonic support is assumed.",
        "canonical_spatial_L2_jet_norm": "U_f(t)=max_{|alpha|<=3} ||partial^alpha f(t,.)||L2(R3), including time and mixed derivatives; integrated norm is ||U_f||L2(time).",
        "fixed_family": "kappa>=kappa0=1e800; fixed COMPLETE canonical functions, mass1000 and zeta1e-6; no bounce assertion away from the anchor",
        "outer_complex_domain": {"abs_u": s.Rational(1, 4), "abs_X": rho},
        "outer_regular_R_minus_one_over_X_squared_bound": 20 * n,
        "inner_regular_coefficient_derivative_bounds": {"u": 160 * n, "X": 640 * n**2},
        "canonical_polydisc": "Each of15 independent Phi, gradient and symmetric Hessian coordinates has modulus<=2; unit-radius polydiscs about real unit jets fit inside.",
        "source_component_holomorphic_upper": M,
        "source_vector_L2_upper_per_U": C,
        "spatial_source_gradient_L2_upper_per_U": 2 * C,
        "source_Frechet_vector_L2_upper_per_U_test": C,
        "spatial_Frechet_gradient_L2_upper_per_U_test": 62 * C,
        "canonical_U_component_terms_upper": (V, B, Ct),
        "checks": checks,
        "gates": {
            "actual_hierarchy_covers_inner_complex_domain": K0 >= 512 * N,
            "switch_order_even_at_least_four": N >= 4 and N % 2 == 0,
            "outer_T_denominator_gap": s.Rational(15, 16) - s.Rational(1, 64) ** 4
            > s.Rational(7, 8),
            "outer_T_over_X_squared_at_most_one": 2 * s.Rational(1, 64) ** 2 < 1,
            "outer_exponential_below_two": 1 / (1 - s.Rational(1, 1024)) < 2,
            "outer_u_denominator_cubed_above_half": s.Rational(15, 16) ** 3
            > s.Rational(1, 2),
            "outer_R_gap": 20 * 4 * s.Rational(1, 64) ** 2 < s.Rational(1, 2),
            "complex_fixed_family_R_gap": s.Rational(5120 * 4, (512 * 4) ** 2)
            < s.Rational(1, 2),
            "uniform_U_below_twenty_thousand": endpoint.subs(n, 4) < 20000,
        },
    }
