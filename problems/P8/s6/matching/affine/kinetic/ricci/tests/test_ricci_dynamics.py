"""Full constraints, coefficient-positive domains and exceptional controls."""
import pytest
import sympy as sp
from p8_affine_ricci import dynamics as d


def test_all_exact_dynamical_identities_and_continuous_sign_proofs():
    for value in d.checks().values():
        assert all(entry == 0 for entry in value) if isinstance(value, sp.MatrixBase) else value == 0
    assert all(value is True for value in d.proof_checks().values())


@pytest.mark.parametrize("time", (sp.Rational(1, 2), -sp.Rational(1, 2), 3))
@pytest.mark.parametrize("coupling", (sp.Rational(1, 1000), sp.Rational(5, 12), 100))
def test_full_actual_lapse_and_shift_solve(time, coupling):
    old = d.old
    bg = old.background()
    values = {symbol: bg[name].subs(d.u, time) for symbol, name in
              ((old.theta, "theta"), (old.ell, "ell"), (old.lam, "lam"), (old.w, "w"), (old.J, "J"))}
    values[old.q] = 17
    values[d.AUXILIARY] = d.coefficients()["A"].subs({d.u: time, d.COUPLING: coupling})
    action = d.scalar_action()["before"].subs(values)
    solution = sp.solve([sp.diff(action, item) for item in (old.n, old.shift)], (old.n, old.shift))
    assert len(solution) == 2
    kinetic = (sp.hessian(action.subs(solution), (old.vd, old.sd))/2).applyfunc(sp.factor)
    assert kinetic == d.scalar_action()["kinetic"].subs(values)
    assert kinetic[1, 1] > 0 and kinetic.det() > 0
    gradient = d.gradients()["matrix"].subs({d.u: time, d.COUPLING: coupling})
    assert gradient[1, 1] > 0 and gradient.det() > 0


def test_strict_gradient_polynomials_cover_the_entire_positive_coupling_domain():
    result = d.gradients()
    ss, determinant = d.positive_polynomial(result["ss_numerator"]), d.positive_polynomial(result["det_numerator"])
    assert (ss["degrees"], ss["monomials"], ss["positive_constant"]) == ((36, 1), 37, 36)
    assert (determinant["degrees"], determinant["monomials"], determinant["positive_constant"]) == ((74, 4), 180, 97119)
    assert ss["strictly_positive"] is determinant["strictly_positive"] is True
    # Gvv itself has negative monomial coefficients: the determinant/ss
    # argument is essential, not a false all-positive-coefficients claim.
    assert d.positive_polynomial(sp.fraction(result["matrix"][0, 0])[0])["strictly_positive"] is False


def test_gyroscopic_term_and_time_dependent_boundary_are_retained():
    result = d.gradients()
    assert result["boundary_residual"] == 0
    assert result["gyroscopic_coefficient"] != 0
    bg = d.old.background()
    wrong_cross = -bg["lam"]*bg["ell"]/(2*bg["theta"])
    difference = sp.factor(result["matrix"][0, 1]-wrong_cross)
    assert difference.subs({d.u: sp.Rational(1, 2), d.COUPLING: 1}) != 0


def test_independent_principal_Hamiltonian_completes_the_gyroscopic_square():
    k00, k01, k11, g00, g01, g11, mixing = sp.symbols("k00 k01 k11 g00 g01 g11 mixing", real=True)
    K = sp.Matrix([[k00, k01], [k01, k11]])
    G = sp.Matrix([[g00, g01], [g01, g11]])
    x, velocity, momentum = (sp.Matrix(sp.symbols(name, real=True)) for name in
                             ("v s", "vd sd", "pv ps"))
    shift = mixing*sp.Matrix([x[1], -x[0]])
    L = (velocity.T*K*velocity)[0]+(velocity.T*shift)[0]-d.q*(x.T*G*x)[0]
    solution = K.inv()*(momentum-shift)/2
    H = (momentum.T*solution)[0]-L.subs(dict(zip(velocity, solution, strict=True)), simultaneous=True)
    target = ((momentum-shift).T*K.inv()*(momentum-shift))[0]/4+d.q*(x.T*G*x)[0]
    assert sp.factor(H-target) == 0


@pytest.mark.parametrize("coupling,time", ((-sp.Rational(1, 1000), 1), (-1, 11), (-3, 19)))
def test_every_negative_coupling_has_a_regular_tail_kinetic_witness(coupling, time):
    witness = d.negative_witness(coupling, time)
    assert witness["denominator"] > sp.Rational(1, 2)
    assert witness["q"] > witness["threshold"] > 0
    assert witness["J_eff"] < 0


@pytest.mark.parametrize("coupling,time", ((-sp.Rational(3, 49), sp.sqrt(sp.Rational(5, 7))),
                                          (-sp.Rational(1, 8), 0), (-sp.Rational(1, 8), sp.sqrt(5))))
def test_singular_curved_algebraic_charts_are_not_inverted(coupling, time):
    with pytest.raises(ValueError, match="nonsingular"):
        d.require_regular(coupling, time, 1, punctured=False)


def test_center_uses_smooth_first_order_coefficients_and_no_Theta_inverse():
    correction = d.coefficients()["A"]
    assert correction.subs(d.u, 0) == 0
    assert sp.diff(correction, d.u).subs(d.u, 0) == 0
    action = d.scalar_action()
    assert not sp.denom(sp.together(action["first_order_H"])).has(d.old.theta)
    assert action["first_order_residual"] == action["first_order_lapse_residual"] == 0
    assert d.require_regular(1, 0, 10, punctured=False)["denominator"] == 9
    with pytest.raises(ValueError, match="u!=0"):
        d.require_regular(1, 0, 10)


@pytest.mark.parametrize("momentum", (sp.Rational(61, 10), 7, 10, 100))
def test_direct_regular_center_kinetic_block(momentum):
    data = d.center()
    kinetic = data["kinetic"].subs(d.q, momentum)
    assert kinetic == data["expected"].subs(d.q, momentum)
    assert kinetic[0, 0] > 0 and kinetic.det() > 0
    assert sp.factor(data["momentum_determinant"]-100*(d.q-6)/(1199*d.q)) == 0


def test_nonunit_normalization_and_missing_commutator_control():
    units = d.units(3, 2, 5)
    assert units["normalized_lambda"] == sp.Rational(5, 12)
    assert units["normalized_C_factor"] == 4
    assert units["additional_mass_or_gap_assigned"] is False
    values = {d.COUPLING: units["normalized_lambda"], d.u: sp.Rational(1, 2)}
    assert d.coefficients()["denominator"].subs(values) == sp.Rational(103, 15)
    wrong = d.COUPLING*d.coefficients()["Cstar_coefficient"]**2/2
    assert sp.factor((wrong-d.coefficients()["A"]).subs(values)) != 0


def test_zero_coupling_and_nonzero_physical_action_mismatch():
    assert d.require_regular(0, 0, 1, punctured=False)["unchanged_auxiliary"] is True
    assert d.coefficients()["A"].subs(d.COUPLING, 0) == 0
    assert d.scalar_action()["physical_quadratic_mismatch"] == 2*d.q*d.AUXILIARY


def test_uniform_momentum_limited_kinetic_comparison_is_not_a_frequency_cutoff():
    result = d.kinetic_remainder()
    assert result["certificate"]["strictly_positive"] is True
    for time in (sp.Rational(1, 100), sp.Rational(1, 2), 3, 10):
        ratio = result["relative_kinetic_norm"].subs({d.u: time, d.COUPLING: sp.Rational(1, 1000), d.q: 2})
        assert 0 < ratio < sp.Rational(1, 50)
    assert "quadratic velocity" in result["scope"]


@pytest.mark.parametrize("bad", (True, False, 1.0, sp.Float(1), "1", sp.I, sp.oo, sp.nan, sp.Symbol("unproved")))
def test_strict_exact_inputs_after_warmup(bad):
    d.require_regular(1, 1, 1)
    with pytest.raises((TypeError, ValueError)):
        d.require_regular(bad, 1, 1)


@pytest.mark.parametrize("values", ((0, 1, 1), (1, 0, 1), (-1, 1, 1), (1, 1.0, 1)))
def test_invalid_physical_units_are_rejected(values):
    with pytest.raises((TypeError, ValueError)):
        d.units(*values)
