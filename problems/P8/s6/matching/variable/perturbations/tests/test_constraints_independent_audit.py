"""Root-authored scalar audit: literal ADM expansion and independent reduction.

These checks do not promote an instantaneous Hamiltonian sign or a
nonuniform frozen momentum limit to a physical health or cone theorem.
"""

from functools import cache

import pytest
import sympy as sp
from p8_variable_constraints import (
    action,
    auxiliary,
    domain,
    independent,
    observables,
    reduction,
    vector,
)


@cache
def literal_diagonal_lagrangian():
    """Expand the lapse/exponential metric densities before any elimination."""
    lam = sp.Symbol("lambda", real=True)
    a, b, c, h, hp, w, kin, b0, b1, b4 = sp.symbols(
        "a b c h hp w kin b0 b1 b4", real=True)
    q, f, e, z = action.QG, action.QF, action.E, action.Z
    vg, vf, ve, vz = action.VG, action.VF, action.VE, action.VZ
    ng, nf, K = action.NG, action.NF, action.K
    sg, sf = 3*q, 3*f+e

    def quadratic(expression):
        return sp.expand(sp.diff(expression, lam, 2).subs(lam, 0)/2)

    # (K_ij K^ij-K^2)/2, including the independent f scalar shear.
    g_eh = -3*a**3*sp.exp(lam*sg)*(h+lam*vg)**2/(1+lam*ng)
    f_eh = (b**3*sp.exp(lam*sf)/(c*(1+lam*nf))
            *(-3*(-h+lam*vf)**2-2*(-h+lam*vf)*lam*ve))
    boundary_g_prime = -a**3*(hp+3*h*h)*sg**2-2*a**3*h*sg*(3*vg)
    boundary_f_prime = b**3*(hp-3*h*h)*sf**2/c+2*b**3*h*sf*(3*vf+ve)/c
    eh_bulk = quadratic(g_eh+f_eh)-boundary_g_prime-boundary_f_prime

    matter = (a**3*sp.exp(lam*sg)/(1+lam*ng)
              *(kin+(w+lam*vz)**2)/2)
    matter_bulk = quadratic(matter)-a*K*z*z/2

    # Literal diagonal principal square root: Nf/Ng, two equal transverse
    # roots and one longitudinal root. beta2=beta3=0 in this action.
    y = b/a
    interaction = -2*a**3*(
        b0*(1+lam*ng)*sp.exp(lam*sg)
        +c*b1*(1+lam*nf)*sp.exp(lam*sg)
        +y*b1*(1+lam*ng)*(2*sp.exp(lam*(2*q+f))+sp.exp(lam*(2*q+f+e)))
        +c*b4*y**3*(1+lam*nf)*sp.exp(lam*sf))
    potential_bulk = quadratic(interaction)
    spatial = a*K*(q*q+2*ng*q)+c*b*K*(f*f+2*nf*f)
    return {"symbols": (a, b, c, h, hp, w, kin, b0, b1, b4),
            "L_no_shift": eh_bulk+matter_bulk+potential_bulk+spatial,
            "EH_boundary_prime": boundary_g_prime+boundary_f_prime}


def test_uneliminated_action_matches_literal_diagonal_adm_density():
    literal, d = literal_diagonal_lagrangian(), action.derive()
    names = ("a", "b", "c", "h", "h_u", "chi_speed", "kbar", "b0", "b1", "b4")
    substitution = dict(zip(literal["symbols"], (d[key] for key in names), strict=True))
    expected = literal["L_no_shift"].subs(substitution, simultaneous=True)
    actual = action.lagrangian()["L"].subs({action.BG: 0, action.BF: 0})
    assert sp.factor(actual-expected) == 0


def test_intrinsic_scalar_gradient_has_no_spurious_longitudinal_shear_term():
    # For gamma_ij=a² exp(2psi) diag(1,1,exp(2E)),
    # sqrt(gamma) R3 = a exp(psi-E)*[-4psi_zz-2psi_z²+4E_z psi_z].
    # Unit Fourier-mode contraction has <cos²>=<sin²>=1. Including lapse
    # and the EH factor 1/2 gives the following quadratic coefficient.
    q, e, ng, a, K = sp.symbols("q e n a K", real=True)
    literal = a*(-4*(-K*q)*(q-e+ng)-2*K*q*q+4*K*e*q)/2
    assert sp.expand(literal-a*K*(q*q+2*ng*q)) == 0


@pytest.mark.parametrize("time,lapse", [
    (0, 4), (0, 3), (0, 1),
    (sp.Rational(1, 100), 4), (sp.Rational(-1, 100), 4),
    (sp.Rational(1, 10), sp.Rational(201, 100)),
])
def test_literal_legendre_and_auxiliary_elimination_gives_H_plus_Nf_Cf(time, lapse):
    d, lag = action.derive(), action.lagrangian()
    at = {d["u"]: time, d["c"]: lapse}
    ll = lag["L"].subs(at)
    velocities = lag["velocities"]
    zero_v = dict.fromkeys(velocities, 0)
    matrix = sp.hessian(ll, velocities)
    linear = sp.Matrix([sp.diff(ll, var).subs(zero_v) for var in velocities])
    momentum = sp.Matrix(d["momenta"])
    hamiltonian = sp.expand(((momentum-linear).T*matrix.inv()*(momentum-linear))[0]/2-ll.subs(zero_v))
    auxiliary = (action.NG, action.BG, action.BF)
    zero_aux = dict.fromkeys(auxiliary, 0)
    aux_matrix = sp.hessian(hamiltonian, auxiliary)
    aux_linear = sp.Matrix([sp.diff(hamiltonian, var).subs(zero_aux) for var in auxiliary])
    reduced = hamiltonian.subs(zero_aux)-(aux_linear.T*aux_matrix.inv()*aux_linear)[0]/2
    assert sp.factor(reduced-(d["H"]+action.NF*d["C_f"]).subs(at)) == 0
    assert sp.factor(sp.diff(reduced, action.NF, 2)) == 0


def test_time_dependent_canonical_map_has_the_stated_one_form_boundary():
    # Independent generic r(u),alpha(u), not a frozen-time symplectic test.
    u = sp.Symbol("u", real=True)
    r, alpha = sp.Function("r")(u), sp.Function("alpha")(u)
    qq, rr, ee, p0, pr, pe = [sp.Function(name)(u) for name in ("Q", "R", "E", "p0", "PR", "PE")]
    old_q = (qq+r*rr+ee/3+alpha*pr, rr+alpha*p0, ee)
    old_p = (p0, pr-r*p0, pe-p0/3)
    old_theta = sum(p*sp.diff(q, u) for p, q in zip(old_p, old_q, strict=True))
    new_theta = p0*sp.diff(qq, u)+pr*sp.diff(rr, u)+pe*sp.diff(ee, u)
    boundary = alpha*p0*pr-r*alpha*p0**2/2
    h_correction = (-sp.diff(r, u)*p0*rr-sp.diff(alpha, u)*p0*pr
                    +(r*sp.diff(alpha, u)-alpha*sp.diff(r, u))*p0**2/2)
    assert sp.expand(old_theta-new_theta-sp.diff(boundary, u)+h_correction) == 0


def test_g_shift_and_actual_clock_observable_cancel_eliminated_momentum():
    a, w, K = sp.symbols("a w K", positive=True)
    z, pe, p0 = sp.symbols("z PE p0", real=True)
    old_pe = pe-p0/3
    p_eg = a**3*w*z-old_pe
    bg = (3*p_eg-p0)/(2*a**3*K)
    assert sp.expand(a*a*bg-3*(a**3*w*z-pe)/(2*a*K)) == 0


def test_regular_primary_constraint_is_not_an_H_quotient():
    d, r = action.derive(), reduction.derive()
    assert sp.factor(r["transformed_C_f"]-r["expected_C_f"]) == 0
    for lapse in (1, 3, 4):
        at = {d["u"]: 0, d["c"]: lapse}
        assert sp.factor(r["expected_secondary_D"].subs(at)) == -sp.Rational(5, 24)
        assert r["alpha"].subs(at) == 0
        assert sp.diff(r["alpha"], d["u"]).subs(at) != 0


def test_secondary_coefficient_from_full_transformed_hessian_is_momentum_independent():
    d, r = action.derive(), reduction.derive()
    direction = r["old_state"].diff(reduction.P0)
    literal = (direction.T*sp.hessian(d["H"], d["state"])*direction)[0]
    literal += sp.diff(r["H_time_correction"], reduction.P0, 2)
    assert sp.factor(literal-r["expected_secondary_D"]) == 0
    assert sp.factor(sp.diff(literal, action.K)) == 0


def test_vector_kinetic_by_relative_and_common_shift_coordinates():
    # Eliminate the common shift first, then the relative shift, rather
    # than solving the production two-equation system at once.
    A, B, C, K = sp.symbols("A B C K", positive=True)
    x, y, common, difference = sp.symbols("x y common difference", real=True)
    lagrangian = A*K*(x-common)**2+B*K*(y-common-difference)**2+C*difference**2
    common_solution = sp.solve(sp.diff(lagrangian, common), common)[0]
    relative_lagrangian = sp.factor(lagrangian.subs(common, common_solution))
    expected_relative = A*B*K*(x-y+difference)**2/(A+B)+C*difference**2
    assert sp.factor(relative_lagrangian-expected_relative) == 0
    difference_solution = sp.solve(sp.diff(relative_lagrangian, difference), difference)[0]
    reduced = sp.factor(relative_lagrangian.subs(difference, difference_solution))
    coefficient = sp.factor(sp.diff(reduced, x, 2)/2)
    primary = vector.schur()
    mapped = primary["kinetic"].subs({primary["A"]: A, primary["B"]: B,
                                        primary["C"]: C, primary["K"]: K})
    assert sp.factor(coefficient-mapped) == 0


def test_vector_canonical_normalization_retains_its_boundary_and_pump():
    u = sp.Symbol("u", real=True)
    kinetic = sp.Function("K")(u)
    freq = sp.Function("frequency_squared")(u)
    canonical = sp.Function("V")(u)
    relative = canonical/sp.sqrt(2*kinetic)
    actual = kinetic*(sp.diff(relative, u)**2-freq*relative**2)
    theta = sp.diff(kinetic, u)/(2*kinetic)
    pump = sp.diff(theta, u)+theta**2
    expected = (sp.diff(canonical, u)**2-(freq-pump)*canonical**2)/2
    boundary = -theta*canonical**2/2
    assert sp.simplify(actual-expected-sp.diff(boundary, u)) == 0


def test_vector_center_negative_control_and_positive_branch_are_not_scalar_health():
    d = vector.derive()
    K = action.K
    at_four = {d["u"]: 0, d["c"]: 4}
    assert sp.factor(d["kinetic"].subs(at_four)-8*K/(3*(K+16))) == 0
    assert sp.factor(d["omega_bare_squared"].subs(at_four)-24-3*K/2) == 0
    assert sp.limit(d["kinetic"].subs({d["u"]: 0, d["c"]: 1}), K, sp.oo) < 0


def test_actual_scalar_observables_are_reconstructed_from_original_metric_and_clock():
    d, r = action.derive(), reduction.derive()
    for time in (0, sp.Rational(1, 100)):
        data = observables.observable_jets(time, 4)
        at = {d["u"]: time, d["c"]: 4}
        old = dict(zip(d["state"], r["old_state"], strict=True))
        xi = (d["a"]**2*d["B_g"]).subs(old, simultaneous=True)
        qf = action.QF.subs(old, simultaneous=True)
        original_observables = sp.Matrix([xi, action.Z+d["chi_speed"]*xi, qf-d["h"]*xi])
        original_observables = original_observables.subs(reduction.Q, 0).subs(at)
        p0 = (data["p0"][0]*sp.Matrix(reduction.STATE))[0]
        rows = original_observables.subs(reduction.P0, p0).jacobian(reduction.STATE)
        assert all(sp.factor(value) == 0 for value in rows-data["O"][0])


def test_center_observable_second_jet_keeps_secondary_momentum_time_derivative():
    d, r = action.derive(), reduction.derive()
    data = observables.observable_jets(0, 4)
    alpha_prime = sp.diff(r["alpha"], d["u"]).subs({d["u"]: 0, d["c"]: 4})
    assert alpha_prime != 0
    expected_third_row = 2*alpha_prime*data["p0"][1]
    assert all(sp.factor(value) == 0 for value in data["O"][2][2:3, :]-expected_third_row)
    assert any(value != 0 for value in expected_third_row)


def test_bounce_finite_K_physical_cauchy_map_is_invertible():
    data = observables.equation(0, 4)
    K = action.K
    expected = 5*(K**2-4*K+192)/(2397*K**2)
    assert sp.factor(data["det_phase_map"]-expected) == 0
    assert sp.expand(K**2-4*K+192-((K-2)**2+188)) == 0
    assert all(sp.factor(value) == 0 for value in data["phase_map"]*data["inverse_phase_map"]-sp.eye(6))


def test_conserved_symplectic_form_is_pulled_back_not_guessed_from_H_momenta():
    for time in (0, sp.Rational(1, 100)):
        data = observables.equation(time, 4)
        inverse = data["inverse_phase_map"].subs(action.K, 1)
        literal = inverse.T*reduction.J*inverse
        assert all(sp.factor(value) == 0 for value in literal-data["symplectic"].subs(action.K, 1))
        assert all(sp.factor(value) == 0 for value in literal+literal.T)


def test_scalar_observable_map_has_a_nonuniform_spatial_frequency_limit():
    center = observables.equation(0, 4)
    nearby = observables.equation(sp.Rational(1, 100), 4)
    K = action.K
    center_num, center_den = sp.fraction(sp.factor(center["det_phase_map"]))
    near_num, near_den = sp.fraction(sp.factor(nearby["det_phase_map"]))
    assert (sp.Poly(center_num, K).degree(), sp.Poly(center_den, K).degree()) == (2, 2)
    assert (sp.Poly(near_num, K).degree(), sp.Poly(near_den, K).degree()) == (3, 2)
    assert center["q_q_bracket"] == sp.zeros(3)
    assert any(value != 0 for value in nearby["q_q_bracket"])
    assert center["symplectic"][3:, 3:] == sp.zeros(3)
    assert any(value != 0 for value in nearby["symplectic"][3:, 3:])


def test_global_commuting_clock_companion_map_has_its_exact_time_boundary():
    u = sp.Symbol("u", real=True)
    t = sp.Function("t")(u)
    xi, z, pxi, pz = [sp.Function(name)(u) for name in ("xi", "Z", "Pxi", "PZ")]
    c0 = sp.Rational(1, 10)
    old_e, old_pe, old_px = t*pxi, c0*z-xi/t, pz+c0*t*pxi
    old_theta = old_pe*sp.diff(old_e, u)+old_px*sp.diff(z, u)
    new_theta = pxi*sp.diff(xi, u)+pz*sp.diff(z, u)
    boundary = c0*t*z*pxi-xi*pxi
    h_time = sp.diff(t, u)*xi*pxi/t
    assert sp.expand(old_theta-new_theta-sp.diff(boundary, u)+h_time) == 0


def test_scalar_shift_terms_from_literal_extrinsic_eigenvalues():
    # Single longitudinal Fourier mode, unit <cos²>=<sin²> contraction.
    # Extrinsic eigenvalues are (v,v,v+w)/N, with
    # v=H+q_t-B_z*q_z and w=e_t-B_z*e_z-B_zz.
    H, K, B, q, e, n, qt = sp.symbols("H K B q e n qt", real=True)
    advection = 6*H*K*B*q+2*H*K*B*e
    velocity_shift = -2*qt*K*B
    volume_lapse = -2*H*K*B*(3*q+e-n)
    assert sp.expand(advection+velocity_shift+volume_lapse+2*K*B*(qt-H*n)) == 0


@pytest.mark.parametrize("time,lapse", [(0, 4), (sp.Rational(1, 100), 4),
                                       (sp.Rational(-1, 100), 3)])
def test_recovered_lapses_and_shifts_satisfy_original_variations(time, lapse):
    d, lag = action.derive(), action.lagrangian()
    recovered = auxiliary.derive(time, lapse)
    at = {d["u"]: time, d["c"]: lapse, action.K: 1}
    state = sp.Matrix(reduction.STATE)
    old = recovered["old_phase"].subs(action.K, 1)*state
    old_prime = recovered["old_phase_prime"].subs(action.K, 1)*state
    values = dict(zip(d["coordinates"], old[:4], strict=True))
    values.update(zip(lag["velocities"], old_prime[:4], strict=True))
    for variable, key in zip(lag["auxiliaries"], ("n_g", "n_f", "B_g", "B_f"), strict=True):
        values[variable] = (recovered[key].subs(action.K, 1)*state)[0]
    ll = lag["L"].subs(at)
    for variable in lag["auxiliaries"]:
        residual = sp.diff(ll, variable).subs(values, simultaneous=True)
        assert sp.factor(residual) == 0
    for index, velocity in enumerate(lag["velocities"]):
        residual = sp.diff(ll, velocity).subs(values, simultaneous=True)-old[4+index]
        assert sp.factor(residual) == 0
    # At H=0 the auxiliary and momentum equations alone cannot see nf.
    # Its recovered value must also satisfy the full old phase evolution.
    j8 = sp.zeros(8)
    j8[:4, 4:] = sp.eye(4)
    j8[4:, :4] = -sp.eye(4)
    hessian = sp.hessian(d["H"], d["state"]).subs(at)
    constraint = sp.Matrix([sp.diff(d["C_f"], variable) for variable in d["state"]]).subs(at)
    exact_rhs = j8*(hessian*old+constraint*values[action.NF])
    assert all(sp.factor(value) == 0 for value in old_prime-exact_rhs)
    omitted_rhs = j8*hessian*old
    assert any(sp.factor(value) != 0 for value in old_prime-omitted_rhs)


def test_secondary_margin_bernstein_cover_has_an_independent_full_replay():
    primary, secondary = domain.derive(), independent.margin_coefficients()
    assert primary["c1_degree"] == (27,)
    assert primary["positive_branch_degree"] == (27, 3)
    assert len(secondary["c1"])+len(secondary["positive"]) == 140
    assert tuple(map(sp.Rational, secondary["c1"])) == primary["c1_coefficients"]
    assert tuple(map(sp.Rational, secondary["positive"])) == primary["positive_branch_coefficients"]
    assert all(value > 0 for group in secondary.values() for value in group)


@pytest.mark.parametrize("bad_order", [1.0, True])
def test_warm_cache_cannot_bypass_exact_jet_order_validation(bad_order):
    reduction.jets(0, 4, 1)
    with pytest.raises((TypeError, ValueError)):
        reduction.jets(0, 4, bad_order)


@pytest.mark.parametrize("momentum_squared", [sp.Integer(1), sp.Rational(7, 2)])
def test_cauchy_determinant_second_jet_by_direct_polynomial_determinant(momentum_squared):
    # Direct determinant of a polynomial matrix, independent of the
    # production Jacobi trace formula for first/second determinant jets.
    data = observables.center_time_map()
    t = sp.Symbol("t", real=True)
    c0, c1, c2 = [matrix.subs(action.K, momentum_squared) for matrix in data["phase_map_jets"]]
    polynomial = sp.Poly(sp.expand((c0+t*c1+t*t*c2/2).det(method="domain-ge")), t)
    assert polynomial.nth(0) == data["det0"].subs(action.K, momentum_squared)
    assert polynomial.nth(1) == 0
    assert 2*polynomial.nth(2) == data["det2"].subs(action.K, momentum_squared)
