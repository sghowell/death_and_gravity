"""Independent full Lorentzian-matrix and variational boundary audits.

The potential oracle uses truncated matrix multiplication and Newton's
definition of every elementary symmetric polynomial. It does not import
the implementation's trace-polynomial expansion or curvature table.
"""
import pytest
import sympy as sp
from p8_general_reduction import basis, bridges, general

ETA = sp.diag(1, -1, -1, -1)


def _convolve(left, right):
    return [sum(left[i]*right[k-i] for i in range(k+1)) for k in range(4)]


def _matrix_convolve(left, right):
    return [sum((left[i]*right[k-i] for i in range(k+1)), sp.zeros(4))
            for k in range(4)]


def _literal_potential(beta, ratio, H):
    root = [ratio*sp.eye(4), ratio*H/2, -ratio*H**2/8, ratio*H**3/16]
    power = [sp.eye(4), sp.zeros(4), sp.zeros(4), sp.zeros(4)]
    traces = [None]
    for _ in range(4):
        power = _matrix_convolve(power, root)
        traces.append([sp.trace(part) for part in power])
    elementary = [[sp.S.One, sp.S.Zero, sp.S.Zero, sp.S.Zero]]
    for n in range(1, 5):
        products = [_convolve(elementary[n-i], traces[i]) for i in range(1, n+1)]
        elementary.append([sum((-1)**(i+1)*products[i-1][k]
                               for i in range(1, n+1))/n for k in range(4)])
    return [-2*sum(beta[n]*elementary[n][k] for n in range(5)) for k in range(4)]


def _root_beta(r, b1, b2, b3):
    return (sp.Rational(7, 5), b1, b2, b3, -(b1+3*r*b2+3*r**2*b3)/r**3)


FIXTURES = (
    (sp.Rational(3, 2), sp.Integer(5), sp.Integer(2), sp.Integer(1)),
    (sp.Integer(2), sp.Integer(3), sp.Integer(-2), sp.Integer(-1)),
    (sp.Integer(1), sp.Integer(1), -sp.Rational(2, 3), sp.Rational(1, 3)),
)
MATRICES = (
    sp.Matrix([[2, 1, -2, 3], [1, -1, 2, -1], [-2, 2, 3, 1], [3, -1, 1, -2]]),
    sp.Matrix([[1, -2, 1, 2], [-2, 4, -1, 3], [1, -1, -3, 2], [2, 3, 2, 1]]),
)


@pytest.mark.parametrize("parameters", FIXTURES)
@pytest.mark.parametrize("covariant", MATRICES)
def test_literal_all_five_beta_full_matrix_expansion(parameters, covariant):
    r, b1, b2, b3 = parameters
    beta = _root_beta(r, b1, b2, b3)
    H = ETA*covariant
    p = b1+2*r*b2+r**2*b3
    z = b1+3*r*b2+2*r**2*b3
    t1, t2, t3 = [sp.trace(H**k) for k in (1, 2, 3)]
    actual = _literal_potential(beta, r, H)
    assert actual[0] == -2*(beta[0]+3*r*b1+3*r**2*b2+r**3*b3)
    assert actual[1] == 0
    assert actual[2] == -r*p*(t2-t1**2)/4
    assert actual[3] == r*(z*t1**3-3*(p+z)*t1*t2+(3*p+2*z)*t3)/24


def test_polarized_full_ten_component_hessian_and_signature_control():
    # Reconstruct the Hessian by independent literal potential evaluations.
    r, b1, b2, b3 = FIXTURES[0]
    beta = _root_beta(r, b1, b2, b3)
    p = b1+2*r*b2+r**2*b3
    entries = [(i, j) for i in range(4) for j in range(i, 4)]
    basis_matrices = []
    for i, j in entries:
        h = sp.zeros(4)
        h[i, j] = h[j, i] = 1
        basis_matrices.append(ETA*h)
    diagonal = [_literal_potential(beta, r, h)[2] for h in basis_matrices]
    hessian = sp.zeros(10)
    for i in range(10):
        hessian[i, i] = 2*diagonal[i]
        for j in range(i):
            hessian[i, j] = hessian[j, i] = (
                _literal_potential(beta, r, basis_matrices[i]+basis_matrices[j])[2]
                - diagonal[i] - diagonal[j])
    assert hessian.det() == 3*(r*p)**10/16
    # All time-space directions are present. A Euclidean-sign replacement
    # would reverse these three diagonal entries and not be this Hessian.
    for pair in ((0, 1), (0, 2), (0, 3)):
        assert hessian[entries.index(pair), entries.index(pair)] == r*p


def test_double_root_has_zero_hessian_but_nonzero_literal_cubic():
    r, b1, b2, b3 = FIXTURES[2]
    beta = _root_beta(r, b1, b2, b3)
    z = sp.Symbol("z", real=True)
    froot = beta[1]+3*z*beta[2]+3*z**2*beta[3]+z**3*beta[4]
    assert sp.expand(froot-(z-1)**2) == 0
    assert sp.diff(froot, z).subs(z, r) == 0
    actual = _literal_potential(beta, r, ETA*MATRICES[0])
    assert actual[1] == actual[2] == 0
    assert actual[3] != 0


@pytest.mark.parametrize("parameters", FIXTURES[:2])
def test_ten_literal_stationary_equations_and_reduced_action(parameters):
    r, b1, b2, b3 = parameters
    p = b1+2*r*b2+r**2*b3
    mf2 = sp.Rational(7, 3)
    B = MATRICES[1]
    entries = [(i, j) for i in range(4) for j in range(i, 4)]
    coordinates = sp.symbols("h0:10")
    h = sp.zeros(4)
    for index, (i, j) in enumerate(entries):
        h[i, j] = h[j, i] = coordinates[index]
    mixed = ETA*h
    density = (-r*p*(sp.trace(mixed**2)-sp.trace(mixed)**2)/4
               + mf2*r**2*sp.trace(ETA*B*mixed)/2)
    solution = mf2*r*(B-ETA*sp.trace(ETA*B)/3)/p
    mapping = {coordinates[index]: solution[i, j] for index, (i, j) in enumerate(entries)}
    assert all(sp.diff(density, coordinate).subs(mapping) == 0 for coordinate in coordinates)
    expected = mf2**2*r**3*(sp.trace((ETA*B)**2)-sp.trace(ETA*B)**2/3)/(4*p)
    assert density.subs(mapping) == expected


@pytest.mark.parametrize("parameters", FIXTURES[:2])
@pytest.mark.parametrize("s,t", ((sp.Rational(2, 3), sp.Rational(-1, 5)),
                                 (sp.S.Zero, sp.Rational(7, 5))))
def test_full_curvature_linear_cubic_against_literal_matrix_oracle(parameters, s, t):
    r, b1, b2, b3 = parameters
    beta = _root_beta(r, b1, b2, b3)
    p = b1+2*r*b2+r**2*b3
    z = b1+3*r*b2+2*r**2*b3
    mf2 = sp.Rational(7, 3)
    alpha = mf2*r/p
    Ric, H = MATRICES
    v = sp.Matrix([3, 1, -1, 1])
    vc = ETA*v
    X = (v.T*vc)[0]
    T = sp.trace(ETA*H)
    L1 = sp.trace((ETA*H)**2)
    V = (vc.T*H*vc)[0]
    R = sp.trace(ETA*Ric)
    G = Ric-ETA*R/2
    GH = sp.trace(ETA*G*ETA*H)
    GH2 = sp.trace(ETA*G*(ETA*H)**2)
    GvHv = (vc.T*G*ETA*H*vc)[0]
    Gvv = (vc.T*G*vc)[0]
    scalar_part = -2*s*H+2*(s**2-t)*v*v.T-s**2*X*ETA

    def cubic(curvature_scale):
        relative = alpha*(curvature_scale*(Ric-ETA*R/6)+scalar_part)
        return _literal_potential(beta, r, ETA*relative)[3]

    # Exact coefficient extraction for a polynomial of degree at most 3.
    literal = (-cubic(2)+8*cubic(1)-8*cubic(-1)+cubic(-2))/12
    table = (
        4*s**2*(3*p+2*z)*GH2
        - 8*s**2*(p+z)*T*GH
        + 4*s*(p*s**2-2*(p+z)*t)*X*GH
        - 8*s*(3*p+2*z)*(s**2-t)*GvHv
        + 8*s*(p+z)*(s**2-t)*T*Gvv
        + 4*(s**2-t)*(2*(p+z)*s**2-p*t)*X*Gvv
        + sp.Rational(4, 3)*(2*p+z)*s**2*R*(L1-T**2)
        - sp.Rational(8, 3)*(2*p+z)*s*(s**2-t)*R*V
        + sp.Rational(4, 3)*s*((p+2*z)*s**2-2*(2*p+z)*t)*R*X*T
        + ((p+2*z)*s**4-4*p*s**2*t)*R*X**2)
    assert literal == r*alpha**3*table/8
    if s == 0:
        assert literal == r*alpha**3*p*t**2*X*Gvv/2


def test_generic_boundary_and_literal_tensor_variation():
    assert len(basis.checks()) == 9
    assert set(basis.checks().values()) == {0}
    test = basis.tensor_boundary()
    assert test["complement_coefficient_at_center"] == 0
    assert test["complement_euler"] != 0
    assert test["full_euler"] == 0
    assert not basis.curvature_boundary()["full_action_changed"]
    assert not basis.curvature_boundary()["Xi_projection_is_full_higher_order_invariant"]


@pytest.mark.parametrize("c", (sp.Rational(201, 100), sp.Rational(5, 2), sp.Integer(4)))
def test_actual_boundary_weight_is_nonzero_after_center_differentiation(c):
    values = basis.actual_center_scaling()
    cs, mf2 = values["symbols"]
    derivative = values["W_theta_center"].subs({cs: c, mf2: 1})
    assert derivative < 0
    assert values["W_center"] == 0
    assert values["normalized_Xi_shift_at_X1"].subs(cs, c) > 0


def test_continuous_old_action_and_actual_clock_bridges():
    assert len(bridges.checks()) == 15
    assert set(bridges.checks().values()) == {0}


def test_deformation_against_full_literal_matrix_potential():
    r, b1, b2, b3 = FIXTURES[0]
    beta = _root_beta(r, b1, b2, b3)
    eta = sp.Rational(7, 11)
    increments = (3*r*eta, -2*eta, eta/r, 0, -eta/r**3)
    changed = tuple(left+right for left, right in zip(beta, increments, strict=True))
    for covariant in MATRICES:
        H = ETA*covariant
        before, after = _literal_potential(beta, r, H), _literal_potential(changed, r, H)
        assert before[:3] == after[:3]
        expected = r*eta*(sp.trace(H)**3-3*sp.trace(H)*sp.trace(H**2)+2*sp.trace(H**3))/24
        assert after[3]-before[3] == expected != 0
    assert general.inverse_at(beta, r)["kappa"] == general.inverse_at(changed, r)["kappa"]


def test_nonboundary_control_by_independent_lapse_scale_variation():
    # g=e^(2 epsilon sigma(t))*diag(1,-e^(2Ht),-e^(2Ht),-e^(2Ht)).
    # Compute all four mixed Schouten eigenvalues from the lapse and
    # logarithmic scale velocities, not the implementation's Weyl rule.
    time, epsilon, hubble = sp.symbols("t epsilon H", real=True)
    sigma = sp.Function("sigma")(time)
    beta_dot = hubble+epsilon*sp.diff(sigma, time)
    beta_ddot = epsilon*sp.diff(sigma, time, 2)
    log_lapse_dot = epsilon*sp.diff(sigma, time)
    inverse_lapse2 = sp.exp(-2*epsilon*sigma)
    p0 = (-2*beta_ddot-beta_dot**2+2*beta_dot*log_lapse_dot)*inverse_lapse2
    pspace = -beta_dot**2*inverse_lapse2
    density = sp.exp(3*hubble*time+4*epsilon*sigma)*(3*p0*pspace**2+pspace**3)
    actual = sp.diff(density, epsilon).subs(epsilon, 0)
    box_boundary = sp.diff(sp.exp(3*hubble*time)*sp.diff(sigma, time), time)
    bulk = 8*hubble**6*sp.exp(3*hubble*time)*sigma
    assert sp.simplify(actual-bulk+6*hubble**4*box_boundary) == 0
    # A nonnegative nonzero smooth compact sigma gives a nonzero bulk
    # integral for H!=0. The exact derivative integrates to zero.
    assert bulk.subs({hubble: 1, sigma: 1}) != 0
    assert general.nonboundary_control()["covariant_metric_Euler_trace"] != 0
