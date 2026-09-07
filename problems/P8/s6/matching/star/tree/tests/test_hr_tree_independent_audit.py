"""Root-authored action/clock/graph audits, independent of the child formulas."""

from functools import cache
from itertools import combinations

import pytest
import sympy as sp


@cache
def literal_edge():
    left = sp.symbols("l0:4", positive=True)
    right = sp.symbols("r0:4", positive=True)
    beta = sp.symbols("b0:5", real=True)
    roots = [r/l for l, r in zip(left, right, strict=True)]
    lag = -2*sp.prod(left)*sum(
        beta[k]*sum(sp.prod(group) for group in combinations(roots, k)) for k in range(5))
    ni, nj, ai, aj = sp.symbols("N_i N_j a_i a_j", positive=True)
    iso = dict(zip(left+right, (ni, ai, ai, ai, nj, aj, aj, aj), strict=True))
    jets = [sp.diff(lag, value).subs(iso) for value in left+right]
    rho_i, rho_j = -jets[0]/ai**3, -jets[4]/aj**3
    ps_i = [value/(ni*ai**2) for value in jets[1:4]]
    ps_j = [value/(nj*aj**2) for value in jets[5:8]]
    return ni, nj, ai, aj, beta, rho_i, rho_j, ps_i, ps_j


def test_tree_literal_eight_coframe_variations_and_two_null_weights():
    ni, nj, ai, aj, beta, rhoi, rhoj, psi, psj = literal_edge()
    y, c = aj/ai, nj/ni
    j = beta[1]+2*beta[2]*y+beta[3]*y**2
    assert all(sp.cancel(p-psi[0]) == 0 for p in psi)
    assert all(sp.cancel(p-psj[0]) == 0 for p in psj)
    assert sp.cancel(rhoi+psi[0]-2*(y-c)*j) == 0
    assert sp.cancel(rhoj+psj[0]-2*(c-y)*j/(c*y**3)) == 0
    assert sp.cancel(ni*ai**3*(rhoi+psi[0])+nj*aj**3*(rhoj+psj[0])) == 0


@cache
def bianchi_pair():
    y, c = sp.symbols("y c", positive=True)
    hi, hj = sp.symbols("H_i H_j", real=True)
    beta = sp.symbols("b0:5", real=True)
    a = beta[0]+3*beta[1]*y+3*beta[2]*y**2+beta[3]*y**3
    b = beta[1]+3*beta[2]*y+3*beta[3]*y**2+beta[4]*y**3
    j = beta[1]+2*beta[2]*y+beta[3]*y**2
    rhoi, rhoj = 2*a, 2*b/y**3
    nulli, nullj = 2*(y-c)*j, 2*(c-y)*j/(c*y**3)
    dyi = y*(c*hj-hi)
    ci = sp.diff(rhoi, y)*dyi+3*hi*nulli
    cj = sp.diff(rhoj, y)*dyi/c+3*hj*nullj
    return y, c, hi, hj, beta, j, ci, cj


def test_tree_bianchi_weight_has_two_lapses_not_one():
    y, c, hi, hj, _, j, ci, cj = bianchi_pair()
    assert sp.factor(ci-6*c*j*(y*hj-hi)) == 0
    assert sp.factor(ci+c**2*y**3*cj) == 0
    assert sp.factor(ci+c*y**3*cj) != 0


def test_tree_flux_is_an_unfactored_velocity_difference():
    y, c, hi, hj, _, j, ci, _ = bianchi_pair()
    ai, ni = sp.symbols("a_i N_i", positive=True)
    aj, nj = ai*y, ni*c
    assert sp.factor(ni**2*ai**3*ci-6*ni*nj*ai**2*j*(aj*hj-ai*hi)) == 0


def test_tree_bianchi_and_null_weight_omissions_are_distinct():
    y, c, hi, hj, beta, _, ci, cj = bianchi_pair()
    point = {y: 2, c: 3, hi: 1, hj: 0, **dict(zip(beta, (0, 1, 0, 0, 0), strict=True))}
    assert ci.subs(point) == -18 and cj.subs(point) == sp.Rational(1, 4)
    assert (ci+c*y**3*cj).subs(point) == -12
    assert 2*(y-c)+c**2*y**3*2*(c-y)/(c*y**3) != 0
    assert (2*(y-c)+c**2*y**3*2*(c-y)/(c*y**3)).subs(point) == 4


def test_tree_zero_einstein_coefficient_still_has_bianchi_identity():
    g, h, hp, rho, rhop, pressure = sp.symbols("G H Hprime rho rhoprime pressure", real=True)
    e0 = 3*g*h**2-rho
    es = g*(2*hp+3*h**2)+pressure
    d_e0 = sp.diff(e0, h)*hp+sp.diff(e0, rho)*rhop
    identity = sp.expand(d_e0+3*h*(e0-es)+rhop+3*h*(rho+pressure))
    assert identity == 0
    assert identity.subs(g, 0) == 0


def incidence(vertex_count, edges):
    matrix = sp.zeros(vertex_count, len(edges))
    for col, (left, right) in enumerate(edges):
        matrix[left, col], matrix[right, col] = 1, -1
    return matrix


def test_tree_oriented_incidence_has_no_flux_kernel():
    for edges in ([(0, 1), (2, 1), (2, 3)], [(2, 0), (2, 1), (3, 2)],
                  [(0, 1), (1, 2), (2, 3), (4, 3)]):
        count = len(edges)+1
        matrix = incidence(count, edges)
        assert matrix.nullspace() == []
        assert sp.ones(1, count)*matrix == sp.zeros(1, count-1)
        for omitted in range(count):
            square = matrix.copy()
            square.row_del(omitted)
            assert abs(square.det()) == 1


def test_cycle_circulation_control_is_not_a_claim_of_a_bounce_solution():
    matrix = incidence(3, [(0, 1), (1, 2), (2, 0)])
    assert matrix*sp.ones(3, 1) == sp.zeros(3, 1)
    assert len(matrix.nullspace()) == 1


def test_tree_fixed_component_full_clock_and_chain_rule():
    gs = sp.symbols("G0:4", nonnegative=True)
    ys = (sp.S.One,)+sp.symbols("y1:4", positive=True)
    cs = (sp.S.One,)+sp.symbols("c1:4", positive=True)
    h, hp = sp.symbols("H Hprime", real=True)
    for selected in ((0,), (0, 1), (0, 1, 3), (0, 1, 2, 3)):
        k = sum(gs[i]*ys[i]**2 for i in selected)
        yp = [y*(c-1)*h for y, c in zip(ys, cs, strict=True)]
        kp = sum(2*gs[i]*ys[i]*yp[i] for i in selected)
        full = 0
        for i in selected:
            # H_i=H/y_i and D_i=(c_i*y_i)^-1 D_root.
            dhi = (hp/ys[i]-h*yp[i]/ys[i]**2)/(cs[i]*ys[i])
            full += -2*cs[i]*ys[i]**4*gs[i]*dhi
        assert sp.factor(full+2*k*hp-h*kp) == 0


def test_tree_nonnegative_kinetic_weights_bound_every_component_coefficient():
    gs = [sp.Rational(2), sp.Rational(0), sp.Rational(7), sp.Rational(3)]
    ys = [sp.Rational(1), sp.Rational(3), sp.Rational(2), sp.Rational(1, 2)]
    rates = [sp.Rational(0), sp.Rational(-5), sp.Rational(2), sp.Rational(-3)]
    for count in range(4):
        for rest in combinations((1, 2, 3), count):
            chosen = (0,)+rest
            k = sum(gs[i]*ys[i]**2 for i in chosen)
            coefficient = sum(gs[i]*ys[i]**2*rates[i] for i in chosen)/k
            assert k >= gs[0] > 0
            assert abs(coefficient) <= max(abs(rate) for rate in rates)


def test_tree_actual_all_vertex_stiff_solution_and_scalar_current():
    t = sp.Symbol("T", positive=True)
    ni, nj, ai, aj, beta, rhoi, rhoj, psi, psj = literal_edge()
    a = t**sp.Rational(1, 3)
    point = {ni: 1, nj: 1, ai: a, aj: a,
             **dict(zip(beta, (-3, 1, 0, 0, -1), strict=True))}
    assert all(sp.simplify(value.subs(point)) == 0 for value in (rhoi, rhoj, *psi, *psj))
    h = sp.diff(a, t)/a
    for g in (0, 1, 2, 7):
        phi = sp.sqrt(sp.Rational(2, 3)*g)*sp.log(t)
        density = sp.diff(phi, t)**2/2
        assert sp.simplify(3*g*h**2-density) == 0
        assert sp.simplify(-2*g*sp.diff(h, t)-2*density) == 0
        assert sp.simplify(sp.diff(a**3*sp.diff(phi, t), t)) == 0


def test_tree_actual_algebraic_attachment_does_not_enter_dynamic_component():
    t = sp.Symbol("T", positive=True)
    ni, nj, ai, aj, beta, rhoi, rhoj, psi, psj = literal_edge()
    a = t**sp.Rational(1, 3)
    # Orient this edge from the de Sitter leaf to the stiff center.
    point = {ni: 1/(3*t), nj: 1, ai: a, aj: a,
             **dict(zip(beta, (0, 1, -sp.Rational(1, 2), 0, sp.Rational(1, 2)), strict=True))}
    assert sp.simplify(rhoi.subs(point)-3) == 0
    assert sp.simplify(psi[0].subs(point)+3) == 0
    assert sp.simplify(rhoj.subs(point)) == sp.simplify(psj[0].subs(point)) == 0
    # Add one beta1-tuned leaf identical to the center, with a copied scalar.
    h = 1/(3*t)
    combined_null = 4/(3*t**2)
    assert sp.simplify(-2*2*sp.diff(h, t)-combined_null) == 0
    assert sp.simplify(2*3*sp.diff(h, t)+combined_null) == -2/(3*t**2)


def test_tree_disconnected_zero_physical_kinetic_control_is_necessary():
    t = sp.Symbol("T", real=True)
    h = sp.diff(sp.log(1+t**2), t)
    assert sp.diff(h, t).subs(t, 0) == 2
    ni, nj, ai, aj, beta, rhoi, rhoj, psi, psj = literal_edge()
    point = {ni: 1, nj: 1, ai: 1+t**2, aj: 1,
             **dict(zip(beta, (0, 0, 0, 0, 0), strict=True))}
    # Every lapse and spatial interaction variation vanishes on these actual fields.
    assert all(sp.simplify(value.subs(point)) == 0 for value in (rhoi, rhoj, *psi, *psj))
    # G_i=0, G_j=1; the j metric is Minkowski and both scalar sources are zero.
    assert 3*0*h**2-rhoi.subs(point) == 0
    assert -2*0*sp.diff(h, t)-3*0*h**2-psi[0].subs(point) == 0
    assert sp.simplify(3*h**2) != 0  # Adding a positive EH coefficient changes it.


def test_tree_positive_part_integrating_factor_identity():
    t, cap = sp.symbols("T C", real=True)
    positive_part = sp.Function("Y")(t)
    assert sp.simplify(sp.diff(sp.exp(-cap*t)*positive_part, t)
                       -sp.exp(-cap*t)*(sp.diff(positive_part, t)-cap*positive_part)) == 0


def test_tree_primary_edge_matches_all_separate_coframe_variations():
    from p8_hr_tree import edge
    ni, nj, ai, aj, beta, rhoi, rhoj, psi, psj = literal_edge()
    betas = (2, -3, 1, 4, -2)
    point = {ni: 2, nj: 5, ai: 3, aj: 7, **dict(zip(beta, betas, strict=True))}
    primary = edge.evaluate(betas, sp.Rational(7, 3), sp.Rational(5, 2))
    expected = {"rho_i": rhoi, "rho_j": rhoj, "pressure_i": psi[0], "pressure_j": psj[0],
                "null_i": rhoi+psi[0], "null_j": rhoj+psj[0]}
    assert all(sp.cancel(value.subs(point)-primary[key]) == 0 for key, value in expected.items())


def test_tree_primary_flux_replay_agrees_with_separate_matrix_solver():
    from p8_hr_tree import graph
    edges = [(0, 1), (2, 1), (2, 3)]
    flux = sp.Matrix([sp.Rational(2, 3), -5, 7])
    divergence = incidence(4, edges)*flux
    assert graph.solve_fluxes(4, edges, list(divergence)) == tuple(flux)
    assert graph.solve_fluxes(1, [], [0]) == ()


def test_tree_zero_weight_intermediate_component_and_clock_normalization():
    from p8_hr_tree import components, graph
    chosen = graph.active_component(4, [(0, 1), (1, 2), (1, 3)], [2, 0, -1])
    assert chosen == (0, 1, 3)
    # Here c_i=N_i/N_root, not the principal speed N_i/(N_root*y_i).
    actual = components.reconstruct([2, 0, 7, 3], [1, 3, 2, "1/2"], [1, 1, 5, "1/3"],
                                    [2, 0, 3, 5], -2, chosen)
    assert actual["A"] == sp.Rational(11, 4)
    assert actual["A_fixed_prime"] == 1
    assert actual["weighted_source"] == sp.Rational(53, 24)
    assert actual["lambda"] == sp.Rational(2, 11)
    assert actual["Hprime"] == -sp.Rational(101, 132)
    assert actual["component_log_rate_bound"] == sp.Rational(4, 3)


@pytest.mark.parametrize("count,edges", [(3, [(0, 1), (1, 2), (2, 0)]),
                                       (2, [(0, 1), (1, 0)]), (2, [(0, 0)]),
                                       (3, [(0, 1)]), (2, [(True, 1)]), (2, [(0, 2)])])
def test_tree_invalid_graph_is_not_silently_treated_as_a_tree(count, edges):
    from p8_hr_tree import graph
    with pytest.raises((ValueError, TypeError)):
        graph.validate_tree(count, edges)


@pytest.mark.parametrize("kwargs", [{"Gs": [-1, 0]}, {"Gs": [0, 1]}, {"component": []},
                                   {"component": [0, 0]}, {"component": [1]},
                                   {"component": [True]}, {"root": True}, {"root": 2},
                                   {"ys": [2, 2]}, {"cs": [2, 3]}, {"nulls": [1, -1]}])
def test_tree_physical_positive_einstein_and_component_domain_guards(kwargs):
    from p8_hr_tree import components
    inputs = {"Gs": [1, 0], "ys": [1, 2], "cs": [1, 3], "nulls": [1, 0],
              "H": 0, "component": [0, 1]}
    with pytest.raises((ValueError, TypeError)):
        components.reconstruct(**(inputs | kwargs))


def test_tree_sharp_endpoint_threshold_boundary_is_not_overclaimed():
    from p8_hr_tree import components
    assert components.cd_endpoint_test("1/2", "8/5")["status"] == "INCONCLUSIVE"
    assert components.cd_endpoint_test("1/2", "3/2")["status"].startswith("EXCLUDED")
