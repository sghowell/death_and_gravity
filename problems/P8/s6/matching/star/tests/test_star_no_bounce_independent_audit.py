"""Independent literal HR-potential, proper-clock and branch-scope audits."""

from fractions import Fraction as Q
from itertools import combinations

import pytest
import sympy as sp
from p8_star import background, branches, potential
from p8_trimetric import model


def zero(value):
    entries = list(value) if isinstance(value, sp.MatrixBase) else [value]
    assert all(sp.simplify(entry) == 0 for entry in entries)


def elementary(values, degree):
    return sum((sp.prod(group) for group in combinations(values, degree)), sp.Integer(0))


def literal_link(n, a, nu, au, betas):
    eigenvalues = [nu/n, au/a, au/a, au/a]
    return -2*n*a**3*sum(b*elementary(eigenvalues, degree) for degree, b in enumerate(betas))


def literal_stresses(betas):
    n, a, nu, au = sp.symbols("n a nu au", positive=True)
    lag = literal_link(n, a, nu, au, betas)
    return {"n": n, "a": a, "nu": nu, "au": au,
            "rho_i": sp.cancel(-sp.diff(lag, n)/a**3),
            "p_i": sp.cancel(sp.diff(lag, a)/(3*n*a**2)),
            "rho_u": sp.cancel(-sp.diff(lag, nu)/au**3),
            "p_u": sp.cancel(sp.diff(lag, au)/(3*nu*au**2))}


def test_four_literal_square_root_eigenvalues_generate_every_HR_coefficient():
    ratio, lapse_ratio = sp.symbols("R N", positive=True)
    b0, b1, b2, b3, b4 = sp.symbols("b0 b1 b2 b3 b4", real=True)
    values = [lapse_ratio, ratio, ratio, ratio]
    polynomial = sum(b*elementary(values, i) for i, b in enumerate((b0, b1, b2, b3, b4)))
    a = b0+3*b1*ratio+3*b2*ratio**2+b3*ratio**3
    b = b1+3*b2*ratio+3*b3*ratio**2+b4*ratio**3
    zero(polynomial-a-lapse_ratio*b)
    zero(sp.diff(a, ratio)-3*(b1+2*b2*ratio+b3*ratio**2))
    zero(b-ratio*sp.diff(b, ratio)/3-(b1+2*b2*ratio+b3*ratio**2))


def test_each_spatial_coframe_variation_before_isotropy_gives_the_same_pressure():
    n, nu = sp.symbols("n nu", positive=True)
    leaves = sp.symbols("a1:4", positive=True)
    roots = sp.symbols("u1:4", positive=True)
    a, au = sp.symbols("a au", positive=True)
    betas = sp.symbols("b0:5", real=True)
    eigenvalues = [nu/n]+[x/y for x, y in zip(roots, leaves, strict=True)]
    lag = -2*n*sp.prod(leaves)*sum(b*elementary(eigenvalues, i) for i, b in enumerate(betas))
    isotropic = dict.fromkeys(leaves, a) | dict.fromkeys(roots, au)
    d = literal_stresses(betas)
    dictionary = {d["n"]: n, d["a"]: a, d["nu"]: nu, d["au"]: au}
    for leaf, root in zip(leaves, roots, strict=True):
        zero(sp.diff(lag, leaf).subs(isotropic)/(n*a**2)-d["p_i"].subs(dictionary))
        zero(sp.diff(lag, root).subs(isotropic)/(nu*au**2)-d["p_u"].subs(dictionary))
    zero(-sp.diff(lag, n).subs(isotropic)/a**3-d["rho_i"].subs(dictionary))
    zero(-sp.diff(lag, nu).subs(isotropic)/au**3-d["rho_u"].subs(dictionary))


def test_two_metric_interaction_null_stresses_have_exact_opposite_weighted_signs():
    betas = sp.symbols("beta0:5", real=True)
    d = literal_stresses(betas)
    r, n = d["au"]/d["a"], d["nu"]/d["n"]
    j = betas[1]+2*betas[2]*r+betas[3]*r**2
    leaf_null = d["rho_i"]+d["p_i"]
    root_null = d["rho_u"]+d["p_u"]
    zero(leaf_null+2*(n-r)*j)
    zero(root_null-2*j*(1-r/n)/r**3)
    zero(leaf_null+n*r**3*root_null)


def test_original_auxiliary_action_is_the_beta3_leaf_subfamily_with_central_endpoint():
    ne, ae, nv, av, nu, au = sp.symbols("ne ae nv av nu au", positive=True)
    pg, pf, bg, bf, big_b = sp.symbols("pg pf bg bf B", real=True)
    e, v, u = sp.diag(ne, ae, ae, ae), sp.diag(nv, av, av, av), sp.diag(nu, au, au, au)
    new = literal_link(ne, ae, nu, au, [bg, 0, 0, pg, 0])
    new += literal_link(nv, av, nu, au, [bf, 0, 0, pf, 0])-2*big_b*nu*au**3
    old = model.potential_density(e, v, u, b=big_b, pg=pg, pf=pf, bg=bg, bf=bf)
    zero(new-old)


def test_Bianchi_factorization_keeps_the_interaction_polynomial_undivided():
    r, c = sp.symbols("R c", positive=True)
    hi, hu = sp.symbols("Hi Hu", real=True)
    b0, b1, b2, b3 = sp.symbols("b0 b1 b2 b3", real=True)
    a = b0+3*b1*r+3*b2*r**2+b3*r**3
    j, n = b1+2*b2*r+b3*r**2, r/c
    density, pressure = 2*a, -2*(a+(n-r)*j)
    di_r = r*(n*hu-hi)
    divergence = sp.diff(density, r)*di_r+3*hi*(density+pressure)
    zero(divergence-6*n*j*(r*hu-hi))
    zero(divergence.subs({b1: 0, b2: 0, b3: 0}))
    zero(divergence.subs({b1: 1, b2: -sp.Rational(1, 2), b3: 0, r: 1}))


def test_dynamic_combination_includes_only_healthy_active_Einstein_coefficients():
    gu = sp.Symbol("Gu", nonnegative=True)
    g, f, r, s, c, d = sp.symbols("G F R S c d", positive=True)
    h, hp, j, ell, nh = sp.symbols("H Hp J L nh", real=True)
    kg, kf = g/r**2, f/s**2
    leaf_g = kg*(hp+(1-c)*h**2)-j*(1-c)/r**3
    leaf_f = kf*(hp+(1-d)*h**2)-ell*(1-d)/s**3
    central = gu*hp+nh/2+j*(1-c)/r**3+ell*(1-d)/s**3
    k = gu+kg+kf
    kp = -2*h*(kg*(1-c)+kf*(1-d))
    zero(2*(leaf_g+leaf_f+central)-(2*k*hp-h*kp+nh))
    assert k.is_positive is True


def test_actual_algebraic_branch_solves_both_Einstein_equations_and_the_canonical_source():
    t = sp.Symbol("T", positive=True)
    d = literal_stresses([0, 1, -sp.Rational(1, 2), 0, sp.Rational(1, 2)])
    a, ni, nu = t**sp.Rational(1, 3), 1/(3*t), sp.Integer(1)
    point = {d["a"]: a, d["au"]: a, d["n"]: ni, d["nu"]: nu}
    rho_i, p_i, rho_u, p_u = [sp.simplify(d[key].subs(point)) for key in ("rho_i", "p_i", "rho_u", "p_u")]
    hi, hu = sp.diff(a, t)/(ni*a), sp.diff(a, t)/a
    field = sp.sqrt(sp.Rational(2, 3))*sp.log(t)
    rho = sp.diff(field, t)**2/2
    zero(3*hi**2-rho_i)
    zero(2*sp.diff(hi, t)/ni+3*hi**2+p_i)
    zero(3*hu**2-rho-rho_u)
    zero(2*sp.diff(hu, t)+3*hu**2+rho+p_u)
    zero(sp.diff(a**3*sp.diff(field, t), t))
    assert (rho_i, p_i, rho_u, p_u) == (3, -3, 0, 0)
    assert hi == 1 and hu == 1/(3*t)


def test_actual_algebraic_branch_falsifies_a_single_all_leg_monotonicity_identity():
    t = sp.Symbol("T", positive=True)
    h, nh = 1/(3*t), 2/(3*t**2)
    correct = 2*sp.diff(h, t)+nh  # Central G_u=1; the leaf has J=0.
    wrong = 2*(1+1)*sp.diff(h, t)+nh
    zero(correct)
    zero(wrong+2/(3*t**2))
    assert wrong.is_negative is True


def test_endpoint_only_auxiliary_center_has_a_genuine_undetermined_bounce_control():
    t = sp.Symbol("T", real=True)
    scale = (1+t**2)**2
    u = sp.diag(1, scale, scale, scale)
    out = model.euler_maps(sp.eye(4), sp.eye(4), u, b=0, pg=0, pf=0, bg=0, bf=0,
                           epsilon=0)
    for key in ("E_e", "E_v", "E_u"):
        zero(out[key])
    h = sp.diff(scale, t)/scale
    assert h.subs(t, 0) == 0 and sp.diff(h, t).subs(t, 0) == 4
    # Adding a positive central EH term would make the zero-source lapse
    # equation 3G_u H²=0 and invalidate this undetermined-geometry control.
    assert (3*h**2).subs(t, 1) == 12


def test_algebraic_leg_positive_lapse_preserves_the_sign_of_its_constant_Hubble_rate():
    t = sp.Symbol("T", real=True)
    c = sp.Function("c", positive=True)(t)
    r = sp.Symbol("R", positive=True)
    hi = sp.Symbol("Hi", real=True)
    h = c*hi/r
    zero(sp.diff(h, t)-sp.diff(c, t)*h/c)
    assert (h/hi).is_positive is True


def test_compact_dynamic_coefficient_bound_is_a_positive_weight_average():
    # K'/2K is -sum K_i(R_i'/R_i)/(G_u+sum K_i). The central
    # constant coefficient contributes zero rate and cannot enlarge it.
    for gu, weights, rates in ((Q(5), [Q(2), Q(3)], [Q(1, 2), Q(-3, 4)]),
                               (Q(0), [Q(2)], [Q(-7, 3)]),
                               (Q(1), [Q(1), Q(4)], [Q(0), Q(2)])):
        k = gu+sum(weights)
        coefficient = -sum(w*v for w, v in zip(weights, rates, strict=True))/k
        cap = max(abs(v) for v in rates)
        assert k > 0 and abs(coefficient) <= cap


def test_positive_part_Gronwall_integrating_factor_has_the_required_sign():
    t = sp.Symbol("T", real=True)
    cap = sp.Symbol("C", nonnegative=True)
    positive = sp.Function("Hplus", nonnegative=True)(t)
    zero(sp.diff(sp.exp(-cap*t)*positive, t)
         -sp.exp(-cap*t)*(sp.diff(positive, t)-cap*positive))


def test_primary_potential_outputs_match_independent_literal_two_metric_variations():
    primary = potential.derive()
    independent = literal_stresses(primary["beta"])
    point = {independent["a"]: 1, independent["au"]: primary["R"],
             independent["n"]: 1, independent["nu"]: primary["N"]}
    for own, new in (("rho_i", "rho_i"), ("p_i", "pressure_i"),
                     ("rho_u", "rho_u"), ("p_u", "pressure_u")):
        zero(independent[own].subs(point)-primary[new])
    r, n = primary["R"], primary["N"]
    zero(primary["bianchi"]-6*n*primary["J"]*(r*primary["H_u"]-primary["H_i"]))


def test_numeric_link_interface_does_not_invert_the_actual_lapse_ratio():
    point = potential.evaluate([2, -3, 1, 2, -1], sp.Rational(3, 2), sp.Rational(2, 5))
    assert point["J"] == sp.Rational(9, 2)
    assert point["rho_i"] == 4 and point["pressure_i"] == sp.Rational(59, 10)
    assert point["rho_u"] == sp.Rational(62, 9) and point["pressure_u"] == -sp.Rational(128, 9)
    assert point["null_i"] == sp.Rational(99, 10) and point["null_u"] == -sp.Rational(22, 3)


def test_exact_polynomial_branch_classification_preserves_positive_double_and_irrational_roots():
    cases = (
        ([0, -2, sp.Rational(3, 2), -1, 0], {(sp.Integer(1), 1), (sp.Integer(2), 1)}),
        ([0, 1, -1, 1, 0], {(sp.Integer(1), 2)}),
        ([0, -2, 0, 1, 0], {(sp.sqrt(2), 1)}),
        ([0, 1, 0, 0, 0], set()),
        ([0, 0, sp.Rational(1, 2), 0, 0], set()),  # R=0 is not regular.
    )
    for betas, roots in cases:
        result = branches.classify(betas)
        assert result["kind"] == "genuine" and not result["identically_zero"]
        assert set(result["positive_roots"]) == roots
    assert branches.classify([7, 0, 0, 0, -9]) == {
        "kind": "endpoint_only", "positive_roots": (), "identically_zero": True}


def test_dynamic_reconstruction_matches_separate_Fraction_clock_and_source_algebra():
    gs, rs, cs = [Q(2), Q(3)], [Q(1), Q(2)], [Q(1, 2), Q(3, 2)]
    h, nh, gu = Q(-2, 3), Q(5, 7), Q(7, 11)
    weights = [g/r**2 for g, r in zip(gs, rs, strict=True)]
    k = gu+sum(weights)
    rp = [r*(1-c)*h for r, c in zip(rs, cs, strict=True)]
    kp = sum(-2*g*v/r**3 for g, r, v in zip(gs, rs, rp, strict=True))
    hp = (h*kp-nh)/(2*k)
    got = background.dynamic_reconstruction(gs, rs, cs, h, nh, central_G=gu)
    assert got["K"] == sp.Rational(k) and got["Kprime"] == sp.Rational(kp)
    assert got["Hprime"] == sp.Rational(hp)
    assert got["Rprimes"] == tuple(sp.Rational(x) for x in rp)
    zero(got["scaled_H_prime"]+sp.Rational(nh)/(2*sp.Rational(k)**sp.Rational(3, 2)))


def test_empty_dynamic_stratum_requires_a_positive_central_Einstein_term():
    with pytest.raises(ValueError):
        background.dynamic_reconstruction([], [], [], 1, 1, central_G=0)
    ordinary = background.dynamic_reconstruction([], [], [], -1, 2, central_G=3)
    assert ordinary["K"] == 3 and ordinary["Kprime"] == 0
    assert ordinary["Hprime"] == -sp.Rational(1, 3) and ordinary["Rprimes"] == ()


@pytest.mark.parametrize("bad", [True, 0.5, sp.Float("0.5"), sp.oo, sp.nan, sp.sqrt(2)])
def test_physical_numeric_interfaces_reject_inexact_nonfinite_and_nonrational_data(bad):
    with pytest.raises((TypeError, ValueError)):
        potential.evaluate([0, 1, 0, 0, 0], bad, 1)
    with pytest.raises((TypeError, ValueError)):
        background.dynamic_reconstruction([bad], [1], [1], 1, 1)


@pytest.mark.parametrize("bad", [0, -1])
def test_principal_clock_and_leaf_Einstein_coefficients_must_be_strictly_positive(bad):
    with pytest.raises(ValueError):
        potential.evaluate([0, 1, 0, 0, 0], 1, bad)
    with pytest.raises(ValueError):
        background.dynamic_reconstruction([1], [1], [bad], 1, 1)
    with pytest.raises(ValueError):
        background.dynamic_reconstruction([bad], [1], [1], 1, 1)


def test_compact_comparison_interface_records_supplied_bounds_without_inventing_them():
    assert background.compact_comparison_bound([Q(1, 7), Q(2, 3)], [Q(5, 4), 0]) == sp.Rational(5, 4)
    assert background.compact_comparison_bound([], []) == 0
    with pytest.raises(ValueError):
        background.compact_comparison_bound([1], [])
    with pytest.raises(ValueError):
        background.compact_comparison_bound([-1], [1])


def test_primary_algebraic_fixture_retains_its_actual_nonzero_false_formula_defect():
    d = branches.algebraic_control()
    t = d["T"]
    assert all(sp.simplify(x) == 0 for x in d["residuals"].values())
    zero(d["false_K_all_defect"]+1/(6*sp.sqrt(2)*t**2))
    assert d["false_K_all_defect"].is_negative is True
    zero(d["H_u"]-1/(3*t))
    zero(d["H_i"]-1)
    zero(d["N"]-3*t)


def test_primary_endpoint_only_control_is_checked_against_actual_full_coframe_equations():
    d = branches.endpoint_control()
    u = sp.diag(d["n_u"], d["a_u"], d["a_u"], d["a_u"])
    out = model.euler_maps(sp.eye(4), sp.eye(4), u, b=0, pg=0, pf=0, epsilon=0)
    for key in ("E_e", "E_v", "E_u"):
        zero(out[key])
    zero(sp.diff(d["a_u"], d["T"])/d["a_u"]-d["H_u"])
    assert d["H_u_prime_at_zero"] == 2
