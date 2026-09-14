"""Whole exact identities and independent nonlinear canonical diagnostics."""

import numpy as np
import pytest
import sympy as s
from p8_vacuum_affine_nonlinear_auxiliary_measure import (
    audit,
    canonical,
    cotangent,
    measure,
)
from scipy.linalg import det
from scipy.optimize import root


@pytest.mark.parametrize("name", list(audit.residuals()))
def test_every_whole_exact_identity(name):
    value = audit.residuals()[name]
    assert all(
        entry == 0
        for entry in (list(value) if isinstance(value, s.MatrixBase) else [value])
    )


@pytest.mark.parametrize(
    "name,call,args", audit.bad_cases(), ids=[row[0] for row in audit.bad_cases()]
)
def test_every_domain_and_scope_rejection(name, call, args):
    with pytest.raises((TypeError, ValueError)):
        call(*args)


def test_whole_counts_fixed_parameters_and_unchanged_frontiers():
    assert (
        len(audit.residuals()),
        audit.scalar_entry_count(),
        len(audit.gates()),
        len(audit.controls()),
        audit.rejected_inputs(),
        len(audit.matching()),
    ) == (53, 451, 37, 9, 294, 113)
    assert all(value is True for value in audit.gates().values())
    assert audit.require_parameters(audit.parameters()) == audit.parameters()
    assert audit.frontier() == audit.previous.frontier()
    assert audit.matching()[:-1] == audit.previous.matching()


@pytest.mark.parametrize("cells", (1, 2, 3, 13, 10**30))
def test_explicit_finite_cell_scope_is_not_a_continuum_limit(cells):
    assert audit.require_finite_cells(cells) == cells
    assert audit.require_branch("one_local_regular_auxiliary_root")


@pytest.mark.parametrize("pivot", (s.Rational(-7, 13), -1, 1, s.Rational(1, 10**400)))
def test_positive_density_does_not_replace_regular_pivot_by_positive_pivot(pivot):
    assert audit.require_pivot(pivot) == pivot


@pytest.mark.parametrize("fourth_root_R", (s.Rational(9, 10), 1, s.Rational(41, 40)))
@pytest.mark.parametrize("lapse", (s.Rational(4, 5), 1, s.Rational(6, 5)))
def test_independent_joint_trace_temporal_legendre_solve(fourth_root_R, lapse):
    # Solve the two original Euler equations together, independently of Kstar/Tstar.
    c = canonical
    values = {
        c.R: fourth_root_R**4,
        c.B: s.Rational(2, 17),
        c.F: -s.Rational(7, 5),
        c.j: s.Rational(1, 101),
        c.N: lapse,
        c.Href: s.Rational(3, 29),
        c.p: s.Rational(5, 23),
        c.G: s.Rational(-2, 31),
        c.pm: s.Rational(2, 7),
        c.ph: s.Rational(-1, 11),
        c.h: s.Rational(1, 13),
        c.n: s.Integer(19),
        c.shear: s.Rational(3, 37),
        c.electric: s.Rational(2, 41),
        c.magnetic: s.Rational(5, 43),
        c.wmass: s.Rational(7, 47),
        c.gm: s.Rational(3, 53),
        c.gh: s.Rational(5, 59),
        c.curvature: s.Rational(-1, 61),
        c.zeta: s.Rational(1, 100),
    }
    lag = c.full_trace()["lagrangian"].subs(values, simultaneous=True)
    equations = (s.diff(lag, c.K) - values[c.p], s.diff(lag, c.T) + values[c.G])
    solution = next(iter(s.linsolve(equations, (c.K, c.T))))
    solved = dict(zip((c.K, c.T), solution))
    direct = lapse * (
        values[c.p] * solution[0]
        - lag.subs(solved, simultaneous=True)
        - values[c.G] * solution[1]
        + c.other.subs(values, simultaneous=True)
    )
    whole = c.full_trace()["whole_reduced_Hamiltonian"].subs(values, simultaneous=True)
    assert s.factor(direct - whole) == 0
    assert all(s.factor(eq.subs(solved, simultaneous=True)) == 0 for eq in equations)
    without_Gauss = direct + lapse * values[c.G] * solution[1]
    assert s.factor(without_Gauss - whole) != 0


@pytest.mark.parametrize("eta_value", (0.05, 0.1, 0.2))
@pytest.mark.parametrize("initial", ((-0.1, -0.2), (-0.2, -0.1), (0, 0)))
def test_independent_nonlinear_root_and_delta_jacobian(eta_value, initial):
    packet = measure.nonlinear_example()
    q, N, T, p, _pN, _pT, lam, eta = packet["finite_example_variables"]
    selected = {
        q: s.Rational(1, 10),
        p: s.Rational(1, 5),
        lam: s.Rational(1, 2),
        eta: eta_value,
    }
    equations = s.lambdify(
        (N, T),
        packet["whole_nonlinear_example_constraints"][2:, :].subs(selected),
        "numpy",
    )
    jacobian = s.lambdify(
        (N, T),
        packet["whole_nonlinear_example_auxiliary_Hessian"].subs(selected),
        "numpy",
    )
    result = root(
        lambda z: np.asarray(equations(*z), dtype=float).ravel(),
        initial,
        jac=lambda z: np.asarray(jacobian(*z), dtype=float),
        tol=1e-12,
    )
    assert result.success and np.linalg.norm(equations(*result.x)) < 1e-13
    D = np.asarray(jacobian(*result.x), dtype=float)
    assert np.min(np.linalg.eigvalsh(D)) > 0.98
    step = 1e-5
    numerical = np.column_stack(
        [
            (
                np.asarray(equations(*(result.x + step * np.eye(2)[i]))).ravel()
                - np.asarray(equations(*(result.x - step * np.eye(2)[i]))).ravel()
            )
            / (2 * step)
            for i in range(2)
        ]
    )
    assert np.linalg.norm(D - numerical) < 1e-10
    density = abs(det(D, check_finite=True))
    assert abs(density / abs(det(numerical, check_finite=True)) - 1) < 1e-10
    assert abs(density - 1) > 0.001


@pytest.mark.parametrize("cells", (1, 2, 3, 4))
def test_whole_finite_cell_nonzero_secondary_bracket_and_dirac_measure(cells):
    q = s.Matrix(s.symbols("q0:" + str(cells), real=True))
    p = s.Matrix(s.symbols("p0:" + str(cells), real=True))
    N = s.Matrix(s.symbols("N0:" + str(cells), real=True))
    T = s.Matrix(s.symbols("T0:" + str(cells), real=True))
    pN = s.Matrix(s.symbols("pN0:" + str(cells), real=True))
    pT = s.Matrix(s.symbols("pT0:" + str(cells), real=True))
    aux = N.col_join(T)
    primary = pN.col_join(pT)
    DD = s.zeros(cells)
    if cells > 2:
        for i in range(cells):
            DD[i, (i + 1) % cells] = 1
            DD[i, (i - 1) % cells] = -1
    Gauss = DD * p
    H = sum(
        (N[i] ** 2 + T[i] ** 2) / 2
        + N[i] * q[i]
        + T[i] * p[i]
        + N[i] * T[i] * q[i] ** 2 / 100
        + N[i] ** 2 * T[i] ** 2 / 400
        - N[i] * T[i] * Gauss[i]
        for i in range(cells)
    )
    # A diagnostic cross-cell potential exercises the arbitrary full D block;
    # it is not asserted to be a spatial-lapse term of the physical action.
    H += sum((N[i] - N[(i + 1) % cells]) ** 2 / 1000 for i in range(cells))
    sample = {
        **{q[i]: s.Rational(i + 1, 100) for i in range(cells)},
        **{p[i]: s.Rational(2 * i - 1, 100) for i in range(cells)},
        **{N[i]: s.Rational(90 + i, 100) for i in range(cells)},
        **{T[i]: s.Rational(i + 1, 100) for i in range(cells)},
    }
    forcing = s.Matrix([s.diff(H, item).subs(sample) for item in aux])
    centered = H - forcing.dot(aux)
    secondary = s.Matrix([s.diff(centered, item) for item in aux])
    assert secondary.subs(sample) == s.zeros(2 * cells, 1)
    Q, P = q.col_join(aux), p.col_join(primary)

    def PB(left, right):
        return sum(
            s.diff(left, x) * s.diff(right, y) - s.diff(left, y) * s.diff(right, x)
            for x, y in zip(Q, P)
        )

    constraints = primary.col_join(secondary)
    bracket = s.Matrix(
        4 * cells,
        4 * cells,
        lambda i, j: PB(constraints[i], constraints[j]).subs(sample),
    )
    D = s.hessian(centered, aux).subs(sample)
    E = bracket[2 * cells :, 2 * cells :]
    zero = s.zeros(2 * cells)
    assert E != zero
    assert bracket == zero.row_join(-D).col_join(D.T.row_join(E))
    assert bracket.det() == D.det() ** 2
    assert D.det() != 0
    inverse = (
        (D.T.inv() * E * D.inv())
        .row_join(D.T.inv())
        .col_join((-D.inv()).row_join(zero))
    )
    assert bracket * inverse == s.eye(4 * cells)
    left = s.Matrix([[PB(q[0], item).subs(sample) for item in constraints]])
    right = s.Matrix([PB(item, p[0]).subs(sample) for item in constraints])
    assert (left * inverse * right)[0] == 0
    assert s.Abs(D.det()) / s.Abs(D.det()) == 1
    if cells > 2:
        assert E[:cells, cells:].is_diagonal() is False


@pytest.mark.parametrize("clock", (-0.2, 0, 0.3))
def test_independent_full_combined_point_boundary_phase_symplecticity(clock):
    packet = cotangent.extended_map()
    boundary = cotangent.primitive_boundary_map()
    q, P = (
        packet["whole_extended_new_position_vector"],
        packet["whole_extended_new_momenta"],
    )
    u = s.Symbol("clock", real=True)
    N = q[13]
    C = s.Function("positive_spatial_conformal_factor", positive=True)(u, N)
    I = s.Function("whole_fixed_primitive", real=True)(u, N)
    replacement = {C: s.exp(u) * (1 + N**2 / 10), I: u * N + N**2 / 5}
    positions = (
        packet["whole_extended_old_position_map"]
        .subs(replacement, simultaneous=True)
        .doit()
    )
    momenta = (
        boundary["whole_combined_inverse_momentum_map"]
        .subs(replacement, simultaneous=True)
        .doit()
    )
    full = s.lambdify((u, *q, *P), positions.col_join(momenta), "numpy")
    point = np.array(
        [
            2,
            0.1,
            0.2,
            1.5,
            0.05,
            1.2,
            0.1,
            0.2,
            -0.1,
            0.03,
            -0.02,
            0.01,
            0.05,
            1.1,
            0.03,
            0.04,
            *np.linspace(-0.2, 0.3, 16),
        ]
    )
    step = 1e-5
    jacobian = np.column_stack(
        [
            (
                np.asarray(full(clock, *(point + step * np.eye(32)[i]))).ravel()
                - np.asarray(full(clock, *(point - step * np.eye(32)[i]))).ravel()
            )
            / (2 * step)
            for i in range(32)
        ]
    )
    omega = np.block(
        [[np.zeros((16, 16)), np.eye(16)], [-np.eye(16), np.zeros((16, 16))]]
    )
    assert np.all(np.isfinite(jacobian))
    assert np.linalg.norm(jacobian.T @ omega @ jacobian - omega) < 1e-7
    # The local NumPy determinant warns even for eye(32); use the independent
    # SciPy routine with explicit finite-input checking for this diagnostic.
    assert abs(det(jacobian, check_finite=True) - 1) < 1e-7
    raw = (
        packet["whole_extended_inverse_momentum_map"]
        .subs(replacement, simultaneous=True)
        .doit()
    )
    difference = s.lambdify((u, *q, *P), momenta - raw, "numpy")
    assert np.linalg.norm(difference(clock, *point)) > 0.1


def test_full_source_hessian_not_a_linear_source_in_new_coordinates():
    packet = measure.source_centering()
    assert packet["whole_auxiliary_source_contact"] != s.zeros(2)
    assert packet["whole_offshell_field_map_second_variation_contact"] != s.zeros(3)
    assert packet["actual_old_reference_classical_constraint_residual"] != 0
    assert "not proved zero or nonzero" in packet["reference_Gaussian_boundary"]
    assert "not" in packet["source_boundary"]


def test_whole_spatial_source_and_profile_dependencies_are_explicit():
    c = canonical
    H = c.full_trace()["whole_reduced_Hamiltonian"]
    for item in (
        c.shear,
        c.electric,
        c.magnetic,
        c.wmass,
        c.gm,
        c.gh,
        c.curvature,
        c.G,
        c.pm,
        c.ph,
        c.j,
        c.n,
    ):
        assert H.has(item) and s.diff(H, item) != 0
    D = c.full_trace()["whole_auxiliary_Hessian"]
    assert all(D.has(item) for item in (s.diff(c.j, c.N), s.diff(c.j, c.N, 2)))
    packet = c.current_data()
    parent = measure.realization.parent
    fixed = parent.quantum.fixed_profile()
    components = (
        parent.profile.PV,
        fixed["p_v_H"] / parent.state.KAPPA,
        fixed["p_v_Phi"] / parent.state.KAPPA,
    )
    assert all(value != 0 for value in components)
    assert s.expand(packet["all_three_fixed_vacuum_constants"] - sum(components)) == 0
    assert packet["actual_initial_full_lapse_Schur_enclosure"][0] > s.Rational(151, 100)


def test_gauge_and_quantum_boundaries_remain_explicit():
    packet = measure.block_measure()
    assert packet["retained_coordinate_count_before_spatial_gauge"] == 11
    assert packet["retained_phase_count_before_spatial_gauge"] == 22
    assert packet["classical_phase_count_after_three_regular_spatial_gauges"] == 16
    assert "unevaluated" in packet["remaining_gauge_factor"]
    assert "Multiple roots" in packet["finite_regulator_reduction"]
    assert "continuum" in packet["finite_regulator_reduction"]
