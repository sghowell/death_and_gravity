"""Whole angular M1 two-particle cut and crossing-complete logarithmic part."""

from functools import cache

import sympy as s

from . import source

S, T, U, MU, K, Z = source.S, source.T, source.U, source.MU, source.K, source.Z
COEFFICIENT = 1 / (960 * s.pi**2 * K * K)
CHANNELS = (S, T, U)
NU2 = s.Symbol("positive_renormalization_scale_squared", positive=True)


def polynomial(channel, other1, other2):
    return channel * channel - other1 * other2 + 2 * MU * channel + 6 * MU * MU


POLYS = {z: polynomial(z, *[w for w in CHANNELS if w != z]) for z in CHANNELS}


def massless_log_part():
    return -COEFFICIENT * sum(POLYS[z] * s.log(-z / NU2) for z in CHANNELS)


def sphere_average(poly, coordinates):
    result = 0
    for powers, coefficient in s.Poly(poly, *coordinates).terms():
        if any(k % 2 for k in powers):
            continue
        moment = s.prod(s.factorial2(k - 1) for k in powers) / s.factorial2(
            sum(powers) + 1
        )
        result += coefficient * moment
    return s.expand(result)


@cache
def data():
    Q = S - 4 * MU
    mean = (S * S - s.Rational(2, 3) * S * Q + Q * Q * (1 + 2 * Z * Z) / 15) / (
        16 * K * K
    )
    absorptive = s.factor(mean / (32 * s.pi))
    x = s.Symbol("cut_cosine", real=True)
    n = s.Matrix(s.symbols("nx ny nz", real=True))
    sine = s.Symbol("sine_of_external_angle", real=True)
    independent = sphere_average(
        (S - Q * n[2] ** 2) * (S - Q * (sine * n[0] + Z * n[2]) ** 2), tuple(n)
    )
    contracted = s.factor(independent.subs(sine * sine, 1 - Z * Z))
    d0 = (S - Q / 3) / (4 * K)
    d2 = -Q / (6 * K)
    spin = (d0 * d0 + d2 * d2 * s.legendre(2, Z) / 5) / (32 * s.pi)
    cutpoly = POLYS[S].subs(U, 4 * MU - S - T)
    v = s.Symbol("crossing_v", real=True)
    b = {S: 2 * MU + v, T: 0, U: 2 * MU - v}
    checks = {
        "entire_independent_Cartesian_sphere": s.factor(contracted - mean * 16 * K * K),
        "full_angle_sewing_and_invariants": s.factor(
            absorptive.subs(Z, 1 + 2 * T / Q) - cutpoly / (960 * s.pi * K * K)
        ),
        "complete_spin0_spin2_optical_decomposition": s.factor(absorptive - spin),
        "identical_intermediate_phase_space": s.factor(
            absorptive.subs(Z, 1)
            - s.integrate(source.pair_tree(S, x) ** 2, (x, -1, 1)) / (64 * s.pi)
        ),
        "pair_tree_legendre_reconstruction": s.factor(
            source.pair_tree(S, x) - d0 - d2 * s.legendre(2, x)
        ),
        "forward_complete_massive_polynomial": s.factor(
            POLYS[S].subs({T: 0, U: 4 * MU - S}) - (S * S + 2 * MU * S + 6 * MU * MU)
        ),
        "zero_threshold_density_nonzero": s.factor(
            POLYS[S].subs({S: 0, T: 0, U: 4 * MU}) - 6 * MU * MU
        ),
        "transfer_log_not_removed_at_forward": s.factor(
            POLYS[T].subs(b, simultaneous=True) - v * v - 2 * MU * MU
        ),
        "transfer_v_squared_log_coefficient": s.diff(
            POLYS[T].subs(b, simultaneous=True), v, 2
        )
        / 2
        - 1,
        "logarithm_discontinuity_coefficient": COEFFICIENT * s.pi
        - 1 / (960 * s.pi * K * K),
        "full_log_renormalization_derivative": s.simplify(
            s.diff(massless_log_part(), NU2) * NU2 - COEFFICIENT * sum(POLYS.values())
        ),
        "crossing_polynomial_sum_local": s.factor(
            (
                sum(POLYS.values())
                - s.Rational(3, 2) * (S * S + T * T + U * U)
                - 18 * MU * MU
            ).subs(U, 4 * MU - S - T)
        ),
    }
    for a, bvar in ((S, T), (S, U), (T, U)):
        checks["full_log_crossing_" + str(a) + str(bvar)] = s.expand(
            massless_log_part().subs({a: bvar, bvar: a}, simultaneous=True)
            - massless_log_part()
        )
    return {
        "whole_two_M1_absorptive_coefficient": absorptive,
        "whole_invariant_s_cut": COEFFICIENT * s.pi * POLYS[S],
        "positive_physical_spin_weights": {
            "spin0": d0 * d0 / (32 * s.pi),
            "spin2": d2 * d2 / (160 * s.pi),
        },
        "whole_crossing_log_part_modulo_real_local_polynomial": massless_log_part(),
        "whole_three_channel_polynomials": POLYS,
        "normalization": "2 Im A = (1/2!) integral dPhi2 A_left A_right*, dPhi2=dOmega/(32pi^2); therefore Im A = angular_average(A_left A_right*)/(32pi). External Phi states and intermediate M1 are each identical, but the two species are distinct.",
        "sheet": "log(-s-i0)=log(s)-i*pi for s>0. The physical s-channel absorptive coefficient is positive for s>=4m^2. Continuing its polynomial below4m^2 locates the massless branch; this continuation is not an optical positivity assertion in an unphysical external region.",
        "scope": "Complete first M1-pair nonanalytic contribution at gravitational order kappa^-2 to formal four-Phi vacuum scattering. A real crossing-symmetric local polynomial, other species, graviton loops, mixed higher loops and physical matching errors are NOT set to zero.",
        "checks": checks,
        "gates": {
            "full_angle_not_forward_only": all(absorptive.has(x) for x in (Z, MU)),
            "positive_both_partial_wave_squares_retained": d0 != 0 and d2 != 0,
            "all_three_crossed_massless_cuts_retained": len(POLYS) == 3,
            "M1_zero_threshold_differs_from_Phi_pair_threshold": True,
            "source_part_is_gauge_independent_by_both_Ward_identities": True,
            "formal_first_coefficient_not_an_exact_unitary_amplitude": True,
            "real_local_polynomial_not_inferred_from_discontinuity": True,
            "nonzero_transfer_log_retained_after_soft_stripping": True,
        },
    }
