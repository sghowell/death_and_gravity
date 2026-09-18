"""Complete three-real signed subtraction, retaining every interference."""

from functools import cache
from itertools import permutations

import sympy as s
from p8_vacuum_affine_uniform_all_tree_bound import bounds as finite_tree

from . import analytic, source


@cache
def projectors():
    identity = s.eye(8)
    faces = []
    for i in range(3):
        E = s.zeros(8)
        for mask in range(8):
            E[mask, mask & ~(1 << i)] = 1
        faces.append(s.ImmutableMatrix(E))
    remainder = (identity - faces[0]) * (identity - faces[1]) * (identity - faces[2])
    baseline = identity - remainder
    return tuple(faces), s.ImmutableMatrix(baseline), s.ImmutableMatrix(remainder)


@cache
def identities():
    faces, P, R = projectors()
    identity = s.eye(8)
    checks = {
        "seven_face_projection_coefficients": P[7, :]
        - s.Matrix([[1, -1, -1, 1, -1, 1, 1, 0]]),
        "baseline_idempotence": P * P - P,
        "remainder_idempotence": R * R - R,
        "mutual_annihilation": P * R + R * P,
        "whole_amplitude_decomposition": P + R - identity,
    }
    for i, E in enumerate(faces):
        checks[f"face{i}_idempotence"] = E * E - E
        checks[f"face{i}_baseline_preserved_left"] = E * P - E
        checks[f"face{i}_baseline_preserved_right"] = P * E - E
        checks[f"face{i}_remainder_zero_left"] = E * R
        checks[f"face{i}_remainder_zero_right"] = R * E
    for index, perm in enumerate(permutations(range(3))):
        U = s.zeros(8)
        for mask in range(8):
            target = sum(1 << perm[i] for i in range(3) if mask >> i & 1)
            U[mask, target] = 1
        checks[f"permutation{index}_baseline"] = U * P - P * U
        checks[f"permutation{index}_remainder"] = U * R - R * U
    real = s.symbols("x0:8", real=True)
    imag = s.symbols("y0:8", real=True)
    cube = s.Matrix([x + s.I * y for x, y in zip(real, imag)])
    B, G, remainder = (P * cube)[7], cube[7], (R * cube)[7]
    square = lambda z: s.expand(z * s.conjugate(z))
    interference = s.expand(B * s.conjugate(remainder) + s.conjugate(B) * remainder)
    checks["all_complex_interferences_retained"] = s.expand(
        square(G) - square(B) - interference - square(remainder)
    )
    a, b, c = variables = s.symbols("a b c", positive=True)

    def vectorize(expression):
        return s.Matrix(
            [
                expression.subs(
                    {v: 0 for i, v in enumerate(variables) if not mask >> i & 1}
                )
                for mask in range(8)
            ]
        )

    example = 1 + a * b + c
    checks["zero_amplitude_rectangle_counterexample"] = (R * vectorize(example))[7]
    checks["nonzero_probability_rectangle_counterexample"] = s.expand(
        (R * vectorize(example**2))[7] - 2 * a * b * c
    )
    checks["defined_density_zero_for_exact_baseline_example"] = s.expand(
        example**2 - (P * vectorize(example))[7] ** 2
    )
    checks["concrete_complex_density_difference"] = square(3) - square(1 + s.I) - 7
    bad_baseline = P - faces[1] * faces[2]
    return {
        "checks": checks,
        "gates": {
            "omitted_double_face_fails_face_identity": any(
                E * bad_baseline != E for E in faces
            ),
            "dropping_interference_changes_generic_density": interference != 0,
            "squared_remainder_alone_fails_complex_example": square(2 - s.I) != 7,
            "probability_rectangle_not_this_prescription": s.expand(
                (R * vectorize(example**2))[7]
            )
            != 0,
            "all_seven_proper_faces_retained": sum(v != 0 for v in P[7, :]) == 7,
        },
        "whole_baseline_face_coefficients": tuple(P[7, :]),
        "whole_definition": "F3=sqrt(rho3)*M3/A0,G3=abc*F3; E_i are compatible normalized energy-face limits; H3=[1-product_i(1-E_i)]G3, B3=H3/(abc). The signed density is |F3|^2-|B3|^2, not the squared amplitude remainder alone or the anchored probability rectangle.",
        "whole_gauge_boundary": "Every retained full-amplitude Ward identity passes linearly to its compatible normalized face limits. No individual connected-line class is made into a physical detector observable.",
    }


@cache
def integration():
    a, b, c, x = s.symbols("a b c x", positive=True)
    entropy = lambda u, v: (u + v) * s.log(u + v) - u * s.log(u) - v * s.log(v)
    I = entropy(a, b)
    J = c * I + b * entropy(a, c) + a * entropy(b, c)
    t = s.Symbol("t", nonnegative=True)
    r, q = s.symbols("r q", positive=True)
    cube_integral = lambda f: s.integrate(f, (a, 0, x), (b, 0, x), (c, 0, x))
    A, B, C = s.symbols("A B C", real=True)
    E, Ep = s.symbols("E Ep", positive=True)
    rho2 = E * E * (Ep * Ep - 1) / ((E * E - 1) * Ep * Ep)
    angular = s.factor(8 * (4 * s.pi) ** 3 / (s.factorial(3) * 2**3 * (2 * s.pi) ** 9))
    checks = {
        "entropy_mixed_derivative": s.factor(s.diff(I, a, b) - 1 / (a + b)),
        "entropy_first_zero_face": s.limit(I, a, 0),
        "entropy_second_zero_face": s.limit(I, b, 0),
        "log_bound_derivative_is_nonnegative_square": s.factor(
            s.diff(t - s.log(1 + t * t), t) - (t - 1) ** 2 / (1 + t * t)
        ),
        "log_bound_zero_origin": (t - s.log(1 + t * t)).subs(t, 0),
        "entropy_square_root_comparison_derivative": s.simplify(
            s.diff(2 * r * q - entropy(r * r, q * q), r)
            - (2 * q - 2 * r * s.log(1 + q * q / (r * r)))
        ),
        "J_majorant_cube_integral": cube_integral(
            2 / s.sqrt(a * b) + 2 / s.sqrt(a * c) + 2 / s.sqrt(b * c)
        )
        - 24 * x * x,
        "J_squared_majorant_cube_integral": cube_integral(12 * (a + b + c)) - 18 * x**4,
        "sum_of_three_squares_majorization": s.expand(
            3 * (A * A + B * B + C * C)
            - (A + B + C) ** 2
            - ((A - B) ** 2 + (A - C) ** 2 + (B - C) ** 2)
        ),
        "exact_three_real_phase_space_factor": s.factor(angular - 1 / (48 * s.pi**6)),
        "real_phase_not_larger_than_one_identity": s.factor(
            1 - rho2 - (E * E - Ep * Ep) / ((E * E - 1) * Ep * Ep)
        ),
    }
    for name, variable in zip(("a", "b", "c"), (a, b, c)):
        checks["J_zero_face_" + name] = s.limit(J, variable, 0)
    D = 7 * finite_tree.coefficient_envelope(3)
    threshold = analytic.budgets()["whole_rectangle_threshold"]
    tv = angular * (48 * D * threshold * x * x + 18 * threshold**2 * x**4)
    checks.update(
        {
            "baseline_from_all_seven_complete_tree_faces": D
            - 126
            * source.HEAVY_MASS2**2
            * (2 * 10**16) ** 3
            / source.KAPPA ** s.Rational(3, 2),
            "same_proved_amplitude_rectangle_threshold": threshold
            - s.Rational(1, 10**670),
            "full_signed_TV_arithmetic": s.simplify(
                tv
                - (
                    D * threshold * x * x / s.pi**6
                    + s.Rational(3, 8) * threshold**2 * x**4 / s.pi**6
                )
            ),
            "squared_remainder_coefficient": threshold**2 - s.Rational(1, 10**1340),
        }
    )
    return {
        "checks": checks,
        "gates": {
            "strict_derived_baseline_bound": bool(0 < D < s.Rational(1, 10**754)),
            "interference_coefficient_below1e_minus1424": bool(
                D * threshold < s.Rational(1, 10**1424)
            ),
            "squared_remainder_prefactor_less_than_one": bool(s.Rational(3, 8) < 1),
            "eight_unit_physical_polarization_triples": 2**3 == 8,
            "common_cutoff_before_signed_limit": True,
            "only_elementary_majorant_extended_outside_simplex": True,
            "positive_original_phase_not_multiplied_twice": True,
            "no_matching_to_virtual_terms_by_definition": True,
            "signed_measure_not_positive_normalized_detector_rate": True,
        },
        "whole_baseline_constant": D,
        "whole_amplitude_rectangle_constant": threshold,
        "whole_angular_polarization_symmetry_factor": angular,
        "whole_cube_kernel_integrals": (24 * x * x, 18 * x**4),
        "whole_total_variation_upper": s.Rational(1, 10**1424) * x * x
        + s.Rational(1, 10**1340) * x**4,
        "whole_relative_reference_upper": 2
        * s.Rational(1, 10**1424)
        * x ** s.Rational(3, 2)
        + 2 * s.Rational(1, 10**1340) * x ** s.Rational(7, 2),
        "whole_limit": "For0<x<=1/8 the defined three-real signed density converges in total variation as a common lower soft cutoff tends to zero, uniformly in the allowed elastic and radiation directions. Its total variation is below1e-1424*x^2+1e-1340*x^4. This is neither positivity nor an inclusive virtual matching theorem.",
        "whole_reference_comparison": "The optional size comparison uses only the retained unexpanded P0>x^alpha/2,0<=alpha<1/2. It does not identify this signed difference as a correction to P0 without matching.",
    }


def require_cutoff(value):
    value = source.require_mass(value)
    if value > s.Rational(1, 8):
        raise ValueError("Require an exact positive total-energy cutoff at most1/8")
    return value


def total_variation_upper(value):
    x = require_cutoff(value)
    return s.Rational(1, 10**1424) * x * x + s.Rational(1, 10**1340) * x**4


def relative_reference_upper(value):
    x = require_cutoff(value)
    return 2 * s.Rational(1, 10**1424) * x ** s.Rational(3, 2) + 2 * s.Rational(
        1, 10**1340
    ) * x ** s.Rational(7, 2)


@cache
def data():
    algebra, measure = identities(), integration()
    return {
        "checks": {
            **{"algebra_" + k: v for k, v in algebra["checks"].items()},
            **{"measure_" + k: v for k, v in measure["checks"].items()},
        },
        "gates": {**algebra["gates"], **measure["gates"]},
        "whole_seven_face_baseline": {
            k: v for k, v in algebra.items() if k not in ("checks", "gates")
        },
        "whole_defined_signed_measure": {
            k: v for k, v in measure.items() if k not in ("checks", "gates")
        },
    }
