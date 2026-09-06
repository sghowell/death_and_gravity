"""Retain the pinned coupled whitening, connection and momentum boundary."""

import sympy as sp
from p8_m1_control import model, oscillator
from p8_m1_loops.background import exact_finite_parameter
from p8_m1_physical import quadratic

from . import reduction as r


def baseline_checks():
    checks = {}
    mapping = {quadratic.H: r.H, quadratic.l: r.L, quadratic.theta: r.THETA,
               quadratic.lam: r.LAMBDA, quadratic.w: r.L*(3*r.DELTA-1), quadratic.J: r.J,
               quadratic.q: r.Q, quadratic.Q[0]: r.QG, quadratic.Q[1]: r.QS,
               quadratic.P[0]: r.PG, quadratic.P[1]: r.PM}
    for chart in ("unitary", "gamma"):
        old = quadratic.symbolic(chart)["density"].subs(mapping, simultaneous=True)
        checks[chart+"_full_old_phase_Hamiltonian"] = sp.factor(old-r.chart_data(chart)["h0"])
    delta_compact = (1-model.x**2)**3/2
    checks["pinned_CD_Lambda_delta"] = sp.expand(model.LAMBDA-(1-3*delta_compact))
    checks["matter_background_density_drift"] = sp.expand(model.derivative(model.l, 1)+3*model.H*model.l)
    checks["fixed_comoving_local_q_drift"] = sp.expand(model.derivative(1/model.z)+6*model.x/model.z)
    return checks


def row_map(coordinate_row, momentum_row, whitening, boundary):
    """Return row r=(rY,rP), with Lensing=a^-3/2*r dot (Y,P).

    This is the full inverse canonical map, not p proportional to Ydot.
    Its time-dependent generator is already in the pinned W and Omega.
    """
    rp = whitening*sp.Matrix(momentum_row)
    ry = whitening.T.inv()*sp.Matrix(coordinate_row)-boundary.T*rp
    return ry, rp


def canonical_generator(potential, connection, lensing_y, lensing_p, q, epsilon):
    """The actual two-scalar first-order representative, without eigenmodes.

    Xdot=A X in fixed time units. P is covariant canonical momentum, not
    Ydot; Omega and the symmetric phase boundary are not interchangeable.
    """
    potential, connection = sp.Matrix(potential), sp.Matrix(connection)
    if potential.shape != (2, 2) or connection.shape != (2, 2):
        raise ValueError("Expected two coupled scalar channels")
    if any(sp.simplify(value) != 0 for value in potential-potential.T):
        raise ValueError("Potential must be symmetric")
    if any(sp.simplify(value) != 0 for value in connection+connection.T):
        raise ValueError("Connection must be antisymmetric")
    row = sp.Matrix.vstack(sp.Matrix(lensing_y), sp.Matrix(lensing_p))
    if row.shape != (4, 1):
        raise ValueError("Expected two coordinate and two momentum row entries")
    symplectic = sp.BlockMatrix([[sp.zeros(2), sp.eye(2)], [-sp.eye(2), sp.zeros(2)]]).as_explicit()
    old = sp.BlockMatrix([[-connection, sp.eye(2)], [-potential, -connection]]).as_explicit()
    correction = -sp.Rational(8, 3)*q**2*symplectic*row*row.T
    return {"old": old, "first_order_correction": correction, "representative": old+epsilon*correction}


def generic_checks():
    t11, t21, t22 = sp.symbols("t11 t21 t22", nonzero=True)
    s11, s12, s22 = sp.symbols("s11 s12 s22", real=True)
    cq = sp.Matrix(sp.symbols("cq1 cq2", real=True))
    cp = sp.Matrix(sp.symbols("cp1 cp2", real=True))
    y, p = sp.Matrix(sp.symbols("Y1 Y2")), sp.Matrix(sp.symbols("P1 P2"))
    a = sp.Symbol("a", positive=True)
    t = sp.Matrix([[t11, 0], [t21, t22]])
    boundary = sp.Matrix([[s11, s12], [s12, s22]])
    ry, rp = row_map(cq, cp, t, boundary)
    old_q = a**(-sp.Rational(3, 2))*t.inv()*y
    old_p = a**(-sp.Rational(3, 2))*t.T*(p-boundary*y)
    direct = (cq.T*old_q+cp.T*old_p)[0]
    normalized = (ry.T*y+rp.T*p)[0]
    inverse_phase = sp.BlockMatrix([[t.inv(), sp.zeros(2)], [-t.T*boundary, t.T]]).as_explicit()
    symplectic = sp.BlockMatrix([[sp.zeros(2), sp.eye(2)], [-sp.eye(2), sp.zeros(2)]]).as_explicit()
    values = {"full_inverse_canonical_lensing_map": sp.factor(a**sp.Rational(3, 2)*direct-normalized),
              "volume_cancels_in_quadratic_Weyl_Hamiltonian": sp.factor(a**3*direct**2-normalized**2)}
    values.update({f"volume_removed_inverse_map_symplectic_{index}": sp.factor(value)
                   for index, value in enumerate(inverse_phase.T*symplectic*inverse_phase-symplectic)})
    return values


def point(chart, x_value, q_fixed, ell_ratio=1):
    """Exact optional point bridge; not used to infer interval bounds.

    Time and momenta are in FIXED centre ell0 units. qbar=s^2*qfixed and
    s=ell/ell0. The scalar gamma chart coordinate has one length weight.
    """
    x_value, q_fixed, ell_ratio = map(exact_finite_parameter, (x_value, q_fixed, ell_ratio))
    if (ell_ratio-sp.Rational(1, 2)).is_nonnegative is not True or (2-ell_ratio).is_nonnegative is not True:
        raise ValueError("Require 1/2<=ell/ell0<=2")
    qbar = ell_ratio**2*q_fixed
    base = oscillator.physical_matrices(chart, x_value, qbar)
    data = oscillator.derive(chart)
    mapping = {model.x: x_value, model.z: 1/qbar, model.l: (1-x_value*x_value)**sp.Rational(11, 2)/10}
    evaluate = lambda value: sp.simplify(value.subs(mapping, simultaneous=True))
    d1, d2, chi = (evaluate(data["factors"][key]) for key in ("d1", "d2", "chi"))
    tbar = sp.Matrix([[sp.sqrt(d1), 0], [chi*sp.sqrt(d2), sp.sqrt(d2)]])
    dimensions = sp.diag(ell_ratio**(-data["factors"]["delta"]), 1)
    t = tbar*dimensions
    root_ratio = sp.sqrt(evaluate(data["factors"]["ratio"]))
    boundary = sp.Matrix([[evaluate(data["momentum_boundary11"]), root_ratio*evaluate(data["momentum_boundary12_factor"])],
                          [root_ratio*evaluate(data["momentum_boundary12_factor"]), evaluate(data["momentum_boundary22"])]])/ell_ratio
    substitutions = {r.H: evaluate(model.H)/ell_ratio, r.THETA: evaluate(model.THETA)/ell_ratio,
                     r.LAMBDA: evaluate(model.LAMBDA), r.DELTA: (1-x_value*x_value)**3/2,
                     r.L: mapping[model.l]/ell_ratio, r.J: evaluate(model.J)/ell_ratio**2, r.Q: q_fixed}
    row = r.chart_data(chart)["lensing_row"].subs(substitutions, simultaneous=True)
    ry, rp = row_map(row[:2, 0], row[2:, 0], t, boundary)
    return {"T_fixed": t, "boundary_fixed": boundary,
            "potential_fixed": base["potential"]/ell_ratio**2,
            "connection_fixed": base["connection"]/ell_ratio,
            "lensing_Y_row": ry.applyfunc(sp.simplify), "lensing_P_row": rp.applyfunc(sp.simplify)}


def controls():
    # An elementary exact map whose omitted symmetric boundary changes H1.
    t = sp.Matrix([[2, 0], [1, 3]])
    boundary = sp.Matrix([[5, 7], [7, 11]])
    ry, rp = row_map([2, 3], [5, 7], t, boundary)
    wrong, _ = row_map([2, 3], [5, 7], t, sp.zeros(2))
    return {"omitting_S_changes_normalized_lensing_squared": sp.expand((ry[0]+rp[1])**2-(wrong[0]+rp[1])**2),
            "fixed_local_momentum_drift_is_not_zero": -6*model.x/model.z,
            "gamma_coordinate_has_nonzero_scale_weight": sp.Integer(1)}
