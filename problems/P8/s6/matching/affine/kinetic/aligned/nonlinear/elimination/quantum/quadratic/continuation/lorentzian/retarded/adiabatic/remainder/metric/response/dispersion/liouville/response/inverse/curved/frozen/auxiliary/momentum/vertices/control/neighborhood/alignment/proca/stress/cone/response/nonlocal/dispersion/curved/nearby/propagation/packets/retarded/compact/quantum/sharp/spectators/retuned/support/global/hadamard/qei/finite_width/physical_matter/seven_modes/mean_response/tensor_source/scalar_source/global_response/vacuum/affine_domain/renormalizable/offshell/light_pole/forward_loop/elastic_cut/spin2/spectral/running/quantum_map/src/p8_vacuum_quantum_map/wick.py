"""Actual cubic derivative-map Wick mixing, with evanescent subtraction kept."""

from functools import cache
from itertools import combinations

import sympy as sp
from p8_offshell_vacuum import jets


def covariance(a, b):
    """Coincident mass-one Gaussian covariance divided by the scalar tadpole."""
    if any(type(v) is not tuple or len(v) != 4 for v in (a, b)):
        raise ValueError("Two native four-dimensional multiindices are required")
    if any(type(n) is not int or n < 0 for v in (a, b) for n in v):
        raise ValueError("Derivative indices must be nonnegative native integers")
    counts = tuple(a[i] + b[i] for i in range(4))
    rank = sum(counts)
    if rank > 4:
        raise ValueError("Only the actual map's rank-four covariance is in scope")
    if any(n % 2 for n in counts):
        return sp.Integer(0)
    pairs = rank // 2
    numerator = sp.Integer((-1) ** (pairs + sum(b)))
    for sign, n in zip(jets.SIGNS, counts):
        numerator *= sp.factorial2(n - 1) * sign ** (n // 2)
    return numerator / sp.prod(4 + 2 * j for j in range(pairs))


def contract(expr):
    """Sum the three pairings of each literal cubic jet monomial."""
    variables = sorted(expr.free_symbols.intersection(jets.BY_SYMBOL), key=str)
    if not variables:
        raise ValueError("A nonzero cubic polynomial in the actual jets is required")
    polynomial = sp.Poly(expr, *variables)
    answer = sp.Integer(0)
    for powers, coefficient in polynomial.terms():
        if sum(powers) != 3:
            raise ValueError("Only the actual field-degree-three Wick map is supported")
        factors = [v for v, power in zip(variables, powers) for _ in range(power)]
        for i, j in combinations(range(3), 2):
            k = 3 - i - j
            answer += (
                coefficient
                * covariance(jets.BY_SYMBOL[factors[i]], jets.BY_SYMBOL[factors[j]])
                * factors[k]
            )
    return sp.expand(answer)


@cache
def data():
    literal = jets.data()
    phi = jets.PHI
    lam, gamma, c = sp.symbols("quartic_lambda quartic_gamma redundant_c", real=True)
    B, d = sp.symbols("free_Box_eigenvalue continued_dimension", real=True)
    epsilon = sp.Symbol("dimensional_epsilon")
    E, X, Z = literal["free_equation"], literal["X"], literal["Z"]
    blocks = {
        "phi_cubed": phi**3,
        "phi_X": phi * X,
        "phi_squared_E": phi**2 * E,
        "Z": Z,
        "phi_K_X": phi * jets.K(X),
        "phi_K_phi_E": phi * jets.K(phi * E),
    }
    expected = {
        "phi_cubed": 3 * phi,
        "phi_X": phi,
        "phi_squared_E": E,
        "Z": jets.box(phi) / 4,
        "phi_K_X": 2 * phi - jets.box(phi),
        "phi_K_phi_E": jets.box(jets.box(phi)) + 2 * jets.box(phi) + phi,
    }
    covariant_blocks = {
        "phi_cubed": sp.Integer(3),
        "phi_X": sp.Integer(1),
        "phi_squared_E": B + 1,
        "Z": B / d,
        "phi_K_X": 2 - 4 * B / d,
        "phi_K_phi_E": (B + 1) ** 2,
    }
    a_d = -c / 2 + lam * (B + 3) + 6 * gamma * B / d - gamma * (B + 1) ** 2 / 2
    assembled = (
        -c * covariant_blocks["phi_cubed"] / 6
        + lam * (2 * covariant_blocks["phi_X"] + covariant_blocks["phi_squared_E"])
        + gamma
        * (
            2 * covariant_blocks["Z"]
            + covariant_blocks["phi_X"]
            + covariant_blocks["phi_cubed"] / 3
        )
        - gamma * covariant_blocks["phi_K_X"]
        - gamma * covariant_blocks["phi_K_phi_E"] / 2
    )
    a4 = sp.expand(a_d.subs(d, 4))
    a_epsilon = sp.diff(a_d.subs(d, 4 - 2 * epsilon), epsilon).subs(epsilon, 0)
    tadpole = sp.exp(sp.EulerGamma * epsilon) * sp.gamma(epsilon - 1) / (16 * sp.pi**2)
    tadpole_series = sp.series(tadpole, epsilon, 0, 1).removeO().expand()
    product_finite = -(a4 + a_epsilon) / (16 * sp.pi**2)
    on_shell = sp.expand(product_finite.subs(B, -1))
    expected_literal = sum(
        a4.coeff(B, n)
        * (phi if n == 0 else jets.box(phi) if n == 1 else jets.box(jets.box(phi)))
        for n in range(3)
    )
    p2 = sp.Symbol("loop_Minkowski_momentum_squared")
    return {
        "covariant_single_contraction_building_blocks": covariant_blocks,
        "dimensionally_continued_map_coefficient": a_d,
        "four_dimensional_coefficient": a4,
        "evanescent_coefficient": a_epsilon,
        "MSbar_tadpole_at_reference_one": tadpole_series,
        "MSbar_finite_linear_composite_mixing": product_finite,
        "on_shell_one_loop_overlap_correction": on_shell,
        "actual_literal_jet_contraction": contract(literal["cubic_field_redefinition"]),
        "literal_jet_count": len(jets.JETS),
        "scheme": "Mass-one light Gaussian propagator; d=4-2epsilon; reference MSbar scale one. Subtract the pole of the product before taking d=4. This is an explicit composite-operator prescription, not a cosmological counterterm transfer.",
        "checks": {
            "covariant_blocks_assemble_actual_R": sp.expand(assembled - a_d),
            **{
                "literal_" + name: sp.expand(contract(expr) - expected[name])
                for name, expr in blocks.items()
            },
            "full_actual_R_Wick_contraction": sp.expand(
                contract(literal["cubic_field_redefinition"]) - expected_literal
            ),
            "actual_covariant_d_four_limit": sp.expand(
                a4
                - (
                    -c / 2
                    + 3 * lam
                    - gamma / 2
                    + (lam + gamma / 2) * B
                    - gamma * B**2 / 2
                )
            ),
            "evanescent_dimension_coefficient": sp.expand(
                a_epsilon - 3 * gamma * B / 4
            ),
            "native_Gamma_tadpole_pole_and_finite_part": sp.expand(
                tadpole_series + (1 / epsilon + 1) / (16 * sp.pi**2)
            ),
            "product_finite_part_not_dimension_four_first": sp.expand(
                product_finite + (a4 + 3 * gamma * B / 4) / (16 * sp.pi**2)
            ),
            "actual_on_shell_overlap": sp.expand(
                on_shell - (c / 2 - 2 * lam + 9 * gamma / 4) / (16 * sp.pi**2)
            ),
            "quadratic_moment_equals_tadpole_mod_scaleless_contact": sp.cancel(
                p2 / (p2 - 1) - 1 - 1 / (p2 - 1)
            ),
            "quartic_moment_equals_tadpole_mod_scaleless_contacts": sp.cancel(
                p2**2 / (p2 - 1) - p2 - 1 - 1 / (p2 - 1)
            ),
            "nonzero_evanescent_on_shell_difference": sp.expand(
                on_shell
                + a4.subs(B, -1) / (16 * sp.pi**2)
                - 3 * gamma / (64 * sp.pi**2)
            ),
        },
    }
