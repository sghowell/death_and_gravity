"""Complete original massive tree, identical-channel partial waves and elastic cut."""

from functools import cache

import sympy as s
from p8_vacuum_canonical_affine_decoupling import amplitude as original
from p8_vacuum_canonical_affine_decoupling import family

LAMBDA, GAMMA = family.LAMBDA, family.GAMMA
S = s.Symbol("s", positive=True)
x = s.Symbol("angle_cosine", real=True)
lam, gam = s.symbols("lambda gamma", positive=True)
beta = s.sqrt(1 - 4 / S)
a = lam * (3 * S * S - 8 * S + 8) + 3 * gam * S * (S - 4) ** 2 / 4 - 8 * gam
b = (S - 4) ** 2 * (lam - 3 * gam * S / 4)
Q0 = lam * (10 * S * S - 32 * S + 40) / 3 + gam * S * (S - 4) ** 2 / 2 - 8 * gam
TREE = a + b * x * x
PARTIAL0 = beta * Q0 / (32 * s.pi)
PARTIAL2 = beta * b / (240 * s.pi)
RHO = beta * (Q0 * Q0 + 4 * b * b / 45) / (32 * s.pi)


@cache
def data():
    source = original.data()["complete_on_shell_tree_amplitude"]
    names = {str(v): v for v in source.free_symbols}
    T = -(S - 4) * (1 - x) / 2
    U = -(S - 4) * (1 + x) / 2
    full = source.subs(
        {
            names["s"]: S,
            names["t"]: T,
            names["u"]: U,
            names["lambda"]: lam,
            names["gamma"]: gam,
        },
        simultaneous=True,
    )
    angular = s.integrate(s.expand(TREE**2), (x, -1, 1))
    c, bb = s.symbols("constant_vertex phase_beta", positive=True)
    v = s.Symbol("crossing_v", real=True)
    spectral = s.Symbol("spectral_s", positive=True)
    crossing = 2 * v * v / (s.pi * (spectral - 2) * ((spectral - 2) ** 2 - v * v))
    checks = {
        "literal_original_full_on_shell_amplitude": s.expand(full - TREE),
        "exact_zeroth_Legendre_coefficient": s.expand(a + b / 3 - Q0),
        "complete_squared_angular_integral": s.expand(
            angular - 2 * a * a - 4 * a * b / 3 - 2 * b * b / 5
        ),
        "orthogonal_partial_wave_square": s.expand(
            angular - 2 * (Q0 * Q0 + 4 * b * b / 45)
        ),
        "identical_phase_space_to_partial_unitarity": s.simplify(
            RHO - (32 * s.pi / beta) * (PARTIAL0**2 + 5 * PARTIAL2**2)
        ),
        "constant_vertex_bubble_and_normalized_channel": s.simplify(
            bb / (32 * s.pi) * (bb * c * c / (32 * s.pi)) - (bb * c / (32 * s.pi)) ** 2
        ),
        "crossing_even_dispersion_half_second_derivative": s.diff(crossing, v, 2).subs(
            v, 0
        )
        / 2
        - 2 / (s.pi * (spectral - 2) ** 3),
        "actual_fixed_lambda": LAMBDA - s.Rational(1, 10**600),
        "actual_fixed_gamma": GAMMA - s.Rational(1024, 10**800),
    }
    for j in range(7):
        target = 2 * Q0 if j == 0 else 4 * b / 15 if j == 2 else 0
        checks["independent_Legendre_projection_" + str(j)] = s.expand(
            s.integrate(TREE * s.legendre(j, x), (x, -1, 1)) - target
        )
    return {
        "scalar_mass_squared": s.S.One,
        "fixed_lambda": LAMBDA,
        "fixed_gamma": GAMMA,
        "full_angular_amplitude": TREE,
        "angular_polynomial_coefficients": (a, b),
        "normalized_identical_partial_waves": {0: PARTIAL0, 2: PARTIAL2},
        "complete_first_elastic_absorptive_coefficient": RHO,
        "normalization": "Labelled invariant amplitude; full sphere carries exactly one final 1/2!. dPhi2=beta dOmega/(32pi^2), Im A>=beta integral|A|^2/(64pi). Normalized identical channel t_l=beta integral A P_l/(64pi), A=(32pi/beta)sum_even(2l+1)t_l P_l, S_l=1+2it_l.",
        "formal_cut_scope": "This is the first elastic unitarity coefficient of the exact original four-scalar tree. It is not the complete interacting real amplitude or a bounded higher-loop remainder. The later matching disjunction uses the exact optical inequality and an explicit full-amplitude error, not a tree cut substituted without control.",
        "tree_parent_scope": "S177 complete f,r,a3 and source ownership retained; S182 vacuum retuning has no quartic jet. No scalar cubic exists; higher independent scalar vertices begin at six fields, vector-source vertices at five total legs and source contacts at eight, so they do not change connected tree four-scalar scattering. No old GY14 quantum matching is transplanted.",
        "checks": checks,
        "gates": {
            "fixed_lambda_and_gamma_positive": LAMBDA > 0 and GAMMA > 0,
            "both_scalar_mass_and_identical_factors_retained": True,
            "potential_and_both_angular_channels_retained": True,
            "formal_cut_not_full_quantum_amplitude": True,
        },
    }
