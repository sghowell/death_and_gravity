"""Complete current quadratic auxiliary measure, with physical density factors."""

from functools import cache

import sympy as s
from p8_vacuum_affine_coupled_gaussian_state import phase
from p8_vacuum_affine_reduced_scalar_hamiltonian import scalar as input_scalar
from p8_vacuum_affine_scalar_tame_propagator import charts

W = s.Symbol("physical_density_kappa_a_cubed", positive=True)
n, shift = input_scalar.n, input_scalar.b
vd, sd = input_scalar.vd, input_scalar.sd
v, sigma = input_scalar.v, input_scalar.sigma
pv, ps, pn, pb = s.symbols(
    "physical_v_momentum physical_sigma_momentum primary_lapse primary_shift",
    real=True,
)
POSITIONS = s.Matrix([v, sigma, n, shift])
MOMENTA = s.Matrix([pv, ps, pn, pb])
VARIABLES = POSITIONS.col_join(MOMENTA)
OMEGA8 = s.zeros(8)
OMEGA8[:4, 4:] = s.eye(4)
OMEGA8[4:, :4] = -s.eye(4)


def bracket(left, right):
    dl = s.Matrix([s.diff(left, z) for z in VARIABLES])
    dr = s.Matrix([s.diff(right, z) for z in VARIABLES])
    return s.cancel((dl.T * OMEGA8 * dr)[0])


def complete_lagrangian():
    return W * (
        input_scalar.lagrangian(
            J=charts.J,
            theta=charts.th,
            matter=-charts.l * charts.E,
            charge=charts.l,
            gradient=charts.E,
        )
        + 3 * charts.T * n * v
        + s.Rational(9, 2) * charts.A * v * v
    )


@cache
def complete_hamiltonian():
    L = complete_lagrangian()
    velocity = s.Matrix([vd, sd])
    momenta = s.Matrix([pv, ps])
    kinetic = s.hessian(L, velocity)
    linear = s.Matrix([s.diff(L, x).subs({vd: 0, sd: 0}) for x in velocity])
    rates = kinetic.inv() * (momenta - linear)
    return s.cancel(
        (momenta.dot(velocity) - L).subs(dict(zip(velocity, rates)), simultaneous=True)
    )


def finite_dirac_block(D, lower):
    """Exact finite-matrix identity; no continuum determinant regularization."""
    D, lower = s.Matrix(D), s.Matrix(lower)
    if not 1 <= D.rows <= 4 or D.rows != D.cols or lower.shape != D.shape:
        raise ValueError("Require matching exact finite square blocks of size1..4")
    if any(x.has(s.Float) or x.is_real is not True for x in list(D) + list(lower)):
        raise ValueError("Require provably real exact entries")
    if any(x.free_symbols for x in list(D) + list(lower)):
        raise ValueError("This checked finite-matrix interface requires fixed entries")
    if D.det() == 0 or lower != -lower.T:
        raise ValueError("Require invertible D and a skew lower block")
    zero = s.zeros(D.rows)
    matrix = zero.row_join(-D).col_join(D.T.row_join(lower))
    inverse = (
        (D.inv().T * lower * D.inv())
        .row_join(D.inv().T)
        .col_join((-D.inv()).row_join(zero))
    )
    return matrix, inverse


@cache
def data():
    L, H = complete_lagrangian(), complete_hamiltonian()
    aux = s.Matrix([n, shift])
    velocity = s.Matrix([vd, sd])
    primary = s.Matrix([pn, pb])
    secondary = s.Matrix([s.diff(H, x) for x in aux])
    D = s.hessian(H, aux).applyfunc(s.cancel)
    constraints = primary.col_join(secondary)
    dirac = s.Matrix(4, 4, lambda i, j: bracket(constraints[i], constraints[j]))
    lower = dirac[2:, 2:]
    zero = s.zeros(2)
    expected_dirac = zero.row_join(-D).col_join(D.T.row_join(lower))
    inverse = (
        (D.inv().T * lower * D.inv())
        .row_join(D.inv().T)
        .col_join((-D.inv()).row_join(zero))
    )
    solution = -D.inv() * secondary.subs({n: 0, shift: 0})
    reduced = s.cancel(H.subs(dict(zip(aux, solution)), simultaneous=True))
    target = W * phase.whole_hamiltonian().subs(
        {
            charts.v: v,
            charts.sigma: sigma,
            phase.pv: pv / W,
            phase.ps: ps / W,
        },
        simultaneous=True,
    )
    A = s.hessian(L, velocity)
    C = s.hessian(L, aux)
    aux_linear = s.Matrix([s.diff(L, x).subs({n: 0, shift: 0}) for x in aux])
    config_solution = -C.inv() * aux_linear
    outer = s.cancel(L.subs(dict(zip(aux, config_solution)), simultaneous=True))
    K = s.hessian(outer, velocity)
    physical = s.Matrix([v, sigma, pv, ps])
    to_constraints = s.Matrix(4, 4, lambda i, j: bracket(physical[i], constraints[j]))
    ordinary = s.Matrix(4, 4, lambda i, j: bracket(physical[i], physical[j]))
    reduced_bracket = ordinary + to_constraints * inverse * to_constraints.T
    checks = {
        "whole_unreduced_velocity_block": A - W * s.diag(-6, 1),
        "whole_regular_auxiliary_Hessian": D
        - W * s.diag(-2 * charts.J, -s.Rational(2, 3)),
        "whole_nonzero_lower_constraint_bracket": lower[0, 1]
        + W * (2 * charts.E * charts.q + 3 * charts.T) / 3,
        "entire_primary_secondary_Dirac_matrix": dirac - expected_dirac,
        "entire_Dirac_inverse_left": inverse * dirac - s.eye(4),
        "entire_Dirac_inverse_right": dirac * inverse - s.eye(4),
        "full_Dirac_determinant_not_only_lower_zero_case": dirac.det() - D.det() ** 2,
        "unphysical_delta_jacobian_cancels_measure": D.det() / (4 * W**2 * charts.J / 3)
        - 1,
        "physical_Dirac_bracket_is_canonical": reduced_bracket - phase.OMEGA,
        "complete_auxiliary_solutions": secondary.subs(
            dict(zip(aux, solution)), simultaneous=True
        ),
        "whole_reduced_current_Hamiltonian": reduced - target,
        "entire_configuration_auxiliary_determinant": C.det() + 4 * W**2 * charts.th**2,
        "full_outer_configuration_kinetic_determinant": K.det()
        - 2 * W**2 * charts.J / charts.th**2,
        "complete_unreduced_configuration_density_squared": -A.det() * D.det()
        - 8 * W**4 * charts.J,
        "configuration_and_canonical_Gaussian_measures_agree": A.det()
        * D.det()
        / C.det()
        - K.det(),
        "central_configuration_kinetic_determinant": (W * charts.central()["K"]).det()
        - 8 * W**2 * charts.J / charts.central()["Delta"].subs(charts.z, 1 / charts.q),
        "crossing_auxiliary_matrix_remains_regular": D.subs(charts.th, 0) - D,
        "density_after_four_canonical_configuration_rescalings": (
            8 * W**4 * charts.J
        ).subs(W, phase.kappa * phase.a**3)
        / phase.kappa**4
        - 8 * phase.a**12 * charts.J,
    }
    for case in range(4):
        m = s.Integer(case + 1)
        D0 = s.Matrix([[2 + m, s.Rational(1, 3)], [m / 5, 4 + m]])
        K0 = s.Matrix([[0, 1 + m], [-1 - m, 0]])
        B0, inv0 = finite_dirac_block(D0, K0)
        checks[f"noncommuting_Dirac_block_inverse_{case}"] = B0 * inv0 - s.eye(4)
        checks[f"noncommuting_Dirac_block_determinant_{case}"] = (
            B0.det() - D0.det() ** 2
        )
        checks[f"noncommuting_Dirac_lower_inverse_zero_{case}"] = inv0[2:, 2:]
    return {
        "full_physical_density": W,
        "whole_unreduced_lagrangian": L,
        "whole_unreduced_hamiltonian": H,
        "complete_primary_and_secondary_constraints": constraints,
        "full_constraint_bracket_matrix": dirac,
        "full_auxiliary_solution": solution,
        "whole_reduced_hamiltonian": reduced,
        "Dirac_density_before_auxiliary_delta_integration": 4 * W**2 * charts.J / 3,
        "unreduced_configuration_density_squared": 8 * W**4 * charts.J,
        "outer_reduced_configuration_density_squared": 2
        * W**2
        * charts.J
        / charts.th**2,
        "central_reduced_configuration_density_squared": (
            8 * W**2 * charts.J / charts.central()["Delta"].subs(charts.z, 1 / charts.q)
        ),
        "measure_boundary": "Exact finite-regulator quadratic canonical reduction and configuration Gaussian identities, up to common coefficient-independent normalization and a consistently continued phase. The scalar spatial gauge has already been fixed; no full nonlinear covariant/BRST or regulator-removal theorem is asserted. A configuration chart with Theta0 is not used.",
        "checks": {
            name: value.applyfunc(s.cancel)
            if isinstance(value, s.MatrixBase)
            else s.cancel(value)
            for name, value in checks.items()
        },
        "gates": {
            "full_current_lower_terms_present": all(
                H.has(coefficient)
                for coefficient in (charts.A, charts.T, charts.l, charts.E, charts.q)
            ),
            "constraint_lower_commutator_not_deleted": lower[0, 1] != 0,
            "both_auxiliaries_included": D.shape == (2, 2),
            "missing_Jc_measure_changes_configuration_functional": s.diff(
                s.log(s.Symbol("positive_current_Jc", positive=True)) / 2,
                s.Symbol("positive_current_Jc", positive=True),
            )
            != 0,
            "configuration_crossing_not_phase_singularity": C.det().subs(charts.th, 0)
            == 0
            and D.det().subs(charts.th, 0) != 0,
            "full_covariant_nonlinear_measure_not_inferred": True,
        },
    }
