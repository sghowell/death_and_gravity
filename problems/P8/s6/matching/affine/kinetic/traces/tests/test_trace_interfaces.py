"""All-equation, literal full-constraint and nonunit interface controls."""
import sympy as sp
from p8_affine_traces import bridges as b
from p8_affine_traces import geometry, rank_two


def test_all_exact_interfaces_and_continuous_premises():
    for value in b.checks().values():
        assert all(item == 0 for item in value) if isinstance(value, sp.MatrixBase) else value == 0
    assert all(value is True for value in b.proof_checks().values())


def test_full_rank_two_auxiliary_complement_not_a_selected_connection_ansatz():
    data = b.two_trace_action()
    assert data["N"].shape == (8, 60)
    assert data["N"].rank() == 8
    assert data["projector"].rank() == 52
    assert data["W"].rank() == 8
    assert data["quadratic_rolling_only"] is True
    # A full basis of the kernel, not an assumed 52-component subset.
    kernel = sp.Matrix.hstack(*data["N"].nullspace())
    assert kernel.shape == (60, 52)
    assert (kernel.T*data["M"]*kernel).det(method="domain-ge") != 0
    assert data["N"]*data["lift"] == sp.eye(8)


def test_generic_full_scalar_embedding_has_no_unchecked_temporal_constraint():
    data = b.scalar_embedding()
    assert data["unchanged_shift"] == 0
    assert data["pure_vector_residual"] == sp.zeros(2)
    assert data["full_auxiliary_rank_residual"] == 0
    assert data["pure_vector_kinetic"] == rank_two.longitudinal()["kinetic"].subs(rank_two.Q, b.old.q)


def test_time_dependent_lapse_reconstruction_and_reduced_Euler_chain_rule():
    t = sp.Symbol("time", real=True)
    v, s, n, shift = [sp.Function(name)(t) for name in ("v", "s", "n", "b")]
    theta, ell = 1+t**2, 1/(2+t**2)
    volume, q = (1+t**2)**3, 2/(1+t**2)**2
    # Independent nonconstant first-derivative action, including lapse dynamics.
    L0 = volume*(sp.diff(v, t)**2/3+sp.diff(s, t)**2/2
                 +q*(sp.diff(n, t)+t*n)**2/2+n*sp.diff(s, t)+n**2*s+v*s)
    full = L0+2*volume*q*shift*(theta*n-sp.diff(v, t)-ell*s/2)

    def euler(L, field, order=1):
        return sum((-1)**k*sp.diff(sp.diff(L, sp.diff(field, t, k)), t, k)
                   for k in range(order+1))

    lapse = (sp.diff(v, t)+ell*s/2)/theta
    En0 = euler(L0, n)
    reconstructed = -En0/(2*volume*q*theta)
    reduced = L0.subs(n, lapse).doit()
    assert sp.factor(euler(full, n).subs(shift, reconstructed).doit()) == 0
    assert sp.factor(euler(full, shift).subs(n, lapse).doit()) == 0
    for field, target in ((v, euler(L0, v)-sp.diff(En0/theta, t)),
                          (s, euler(L0, s)+ell*En0/(2*theta))):
        difference = euler(reduced, field, 2)-target.subs(n, lapse).doit()
        assert sp.factor(sp.together(difference)) == 0


def test_nonunit_interface_and_actual_on_clock_scope():
    result = b.calibration()
    assert result["nonunit"]["normalized_zeta"] == sp.Rational(5, 12)
    assert result["nonunit"]["isolated_mass_squared_normalized"] == sp.Rational(32, 5)
    assert result["regular_complement_claim_is_quadratic_rolling_only"] is True
    assert geometry.schur()["D_clock"] == sp.kronecker_product(rank_two.D, geometry.ETA)
