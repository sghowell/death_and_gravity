"""Root-authored determinant, all-coframe, physical-clock and domain checks."""

from functools import cache

import pytest
import sympy as sp


@cache
def literal_coframes():
    fields = [sp.symbols(f"e{i}_0:4", positive=True) for i in range(3)]
    beta = sp.symbols("b0:3", real=True)
    lam = sp.Symbol("lambda", positive=True)
    sums = [sum(beta[i]*fields[i][j] for i in range(3)) for j in range(4)]
    lag = -lam*sp.prod(sums)
    aa, nn = sp.symbols("a0:3", positive=True), sp.symbols("N0:3", positive=True)
    iso = {fields[i][j]: nn[i] if j == 0 else aa[i] for i in range(3) for j in range(4)}
    rho, pressures = [], []
    for i in range(3):
        rho.append(-sp.diff(lag, fields[i][0]).subs(iso)/aa[i]**3)
        pressures.append([sp.diff(lag, fields[i][j]).subs(iso)/(nn[i]*aa[i]**2) for j in (1, 2, 3)])
    return lam, beta, aa, nn, rho, pressures


def test_determinant_twelve_independent_coframe_variations():
    lam, beta, aa, nn, rho, ps = literal_coframes()
    sa, sn = sum(b*a for b, a in zip(beta, aa, strict=True)), sum(b*n for b, n in zip(beta, nn, strict=True))
    for i in range(3):
        assert sp.expand(rho[i]-lam*beta[i]*sa**3/aa[i]**3) == 0
        assert all(sp.expand(p+lam*beta[i]*sn*sa**2/(nn[i]*aa[i]**2)) == 0 for p in ps[i])
    assert sp.factor(sum(nn[i]*aa[i]**3*(rho[i]+ps[i][0]) for i in range(3))) == 0


def test_determinant_unfactored_bianchi_and_off_shell_noether_weights():
    lam, beta, aa, nn, rho, ps = literal_coframes()
    vv = sp.symbols("v0:3", real=True)
    sa, sn, sadot = (sum(b*x for b, x in zip(beta, values, strict=True)) for values in (aa, nn, vv))
    balances = []
    for i in range(3):
        balance = sum(sp.diff(rho[i], aa[j])*vv[j] for j in range(3))/nn[i]
        balance += 3*vv[i]/(nn[i]*aa[i])*(rho[i]+ps[i][0])
        expected = 3*lam*beta[i]*sa**2/(nn[i]*aa[i]**3)*(sadot-sn*vv[i]/nn[i])
        assert sp.factor(balance-expected) == 0
        balances.append(balance)
    assert sp.factor(sum(nn[i]**2*aa[i]**3*balances[i] for i in range(3))) == 0


def test_determinant_full_matrix_auxiliary_stationary_map_sixteen_directions():
    # No diagonal-matrix assumption is made for the auxiliary variation.
    total = sp.Matrix([[2, sp.Rational(1, 3), 0, 0], [0, 3, sp.Rational(1, 5), 0],
                       [0, 0, 4, sp.Rational(1, 7)], [sp.Rational(1, 11), 0, 0, 5]])
    z = sp.Symbol("z", real=True)
    for row in range(4):
        for col in range(4):
            variation = sp.zeros(4)
            variation[row, col] = 1
            w = total+z*variation
            # lambda=1, B=-3/2, Q=U/2; use adjugate, not an assumed inverse rule.
            action = 3*w.det()-sp.trace(w.adjugate()*total)
            assert sp.diff(action, z).subs(z, 0) == 0
            assert action.subs(z, 0) == -total.det()


def test_determinant_fixed_physical_clock_combination():
    g = sp.symbols("G0:3", positive=True)
    y = (sp.S.One,)+sp.symbols("y1:3", positive=True)
    c = (sp.S.One,)+sp.symbols("c1:3", positive=True)
    h, hp = sp.symbols("H Hprime", real=True)
    k = sum(gi*yi**2 for gi, yi in zip(g, y, strict=True))
    yp = [yi*(ci-1)*h for yi, ci in zip(y, c, strict=True)]
    kp = sum(2*g[i]*y[i]*yp[i] for i in range(3))
    geometric = sum(-2*g[i]*c[i]*y[i]**4*(hp/y[i]-h*yp[i]/y[i]**2)/(c[i]*y[i]) for i in range(3))
    assert sp.factor(geometric+2*k*hp-h*kp) == 0


def test_determinant_actual_mixed_sign_stiff_source_full_equations():
    t = sp.Symbol("T", positive=True)
    lam, beta, aa, nn, rho, ps = literal_coframes()
    ys, gs, cc = (1, 2, 3), (1, 2, 3), (-8, sp.Rational(1, 2), -sp.Rational(8, 81))
    scales = [y*t**sp.Rational(1, 3) for y in ys]
    point = {lam: 1, **dict(zip(beta, (1, -1, 1), strict=True)),
             **dict(zip(aa, scales, strict=True)), **dict(zip(nn, ys, strict=True))}
    for i in range(3):
        h = sp.diff(scales[i], t)/(ys[i]*scales[i])
        phi = sp.sqrt(sp.Rational(2*gs[i], 3))*sp.log(t)
        matter_rho = sp.diff(phi, t)**2/(2*ys[i]**2)
        assert sp.simplify(gs[i]*(3*h**2-cc[i])-matter_rho-rho[i].subs(point)) == 0
        assert sp.simplify(gs[i]*(-2*sp.diff(h, t)/ys[i]-3*h**2+cc[i])-matter_rho-ps[i][0].subs(point)) == 0
        assert sp.simplify(sp.diff(scales[i]**3*sp.diff(phi, t)/ys[i], t)) == 0
    assert sum(g*y**2 for g, y in zip(gs, ys, strict=True)) == 36


def test_determinant_disconnected_field_invalidates_wrong_all_metric_K():
    # beta=(1,0): genuine exact independent de Sitter metrics, not a kinematic jet.
    t = sp.Symbol("T", real=True)
    a0, a1 = sp.exp(t), sp.exp(-t)
    h0, h1 = sp.diff(a0, t)/a0, sp.diff(a1, t)/a1
    # lambda=G0=G1=1, Lambda0=2, Lambda1=3, no actual matter.
    assert 3*h0**2-2-1 == 0 and -3*h0**2+2+1 == 0
    assert 3*h1**2-3 == 0 and -3*h1**2+3 == 0
    correct_k, wrong_k = 1, 1+(a1/a0)**2
    assert sp.diff(h0/sp.sqrt(correct_k), t) == 0
    assert sp.diff(h0/sp.sqrt(wrong_k), t).subs(t, 0) == sp.sqrt(2)/2


def test_determinant_singular_lapse_bianchi_control_is_not_a_solution():
    from p8_determinant import controls, potential
    control = controls.singular_sum_bianchi_control()
    assert (control["S_a"], control["S_N"], control["S_a_dot"]) == (2, 0, 0)
    assert control["bianchi"] == (0, 0, 0) and len(set(control["Q_i"])) == 3
    assert control["actual_solution"] is control["branch_health_claim"] is False
    values = potential.evaluate((1, 1, -1), (1, 2, 1), (1, 1, 2))
    assert values["regular_sum"] is False
    assert values["rho"] == (8, 1, -8) and values["pressure"] == (0, 0, 0)
    # G_i=1, no cosmological terms or matter: the displayed velocities fail lapse equations.
    hs = (1, -sp.Rational(1, 2), 0)
    assert tuple(3*h**2-r for h, r in zip(hs, values["rho"], strict=True)) == (-5, -sp.Rational(1, 4), 8)


def test_determinant_primary_literal_and_clock_bridges():
    from p8_determinant import background, potential
    lam, beta, aa, nn, rho, ps = literal_coframes()
    betas, scales, lapses = (2, -1, 3), (2, 3, 5), (3, 2, 7)
    point = {lam: 2, **dict(zip(beta, betas, strict=True)),
             **dict(zip(aa, scales, strict=True)), **dict(zip(nn, lapses, strict=True))}
    primary = potential.evaluate(betas, scales, lapses, lambda_=2)
    assert tuple(sp.factor(r.subs(point)) for r in rho) == primary["rho"]
    assert tuple(sp.factor(p[0].subs(point)) for p in ps) == primary["pressure"]
    rates = background.reconstruct([1, 2, 3], [1, 2, 3], [1, 1, 1], 0, [2, 3, 5])
    assert rates["K"] == 36 and rates["Kprime"] == 0
    assert rates["weighted_null"] == 455
    assert rates["Hprime"] == -sp.Rational(455, 72)


def test_determinant_actual_controls_and_source_preserving_map():
    from p8_determinant import controls
    actual = controls.proportional_stiff()
    assert actual["cosmological_constants"] == (-8, sp.Rational(1, 2), -sp.Rational(8, 81))
    assert actual["K"] == 36
    assert all(value == 0 for value in actual["residuals"].values())
    assert all(value == 0 for value in controls.proportional_de_sitter()["residuals"].values())
    mapping = controls.auxiliary_map()
    assert mapping["matter_stays_on_original_EH_leaves"] is True
    assert mapping["matter_on_auxiliary_w"] is False
    assert all(value == 0 for value in mapping["stationary_Euler"])


@pytest.mark.parametrize("value", [True, 0.1, sp.Float("0.1"), sp.oo, sp.nan, sp.Symbol("x")])
def test_determinant_public_exact_domain_guards(value):
    from p8_determinant import potential
    with pytest.raises((TypeError, ValueError)):
        potential.evaluate((1, value), (1, 2), (1, 2))


@pytest.mark.parametrize("kwargs", [
    {"Gs": [0]}, {"Gs": [-1]}, {"ys": [2]}, {"cs": [2]},
    {"nulls": [-1]}, {"physical": True}, {"physical": 1},
])
def test_determinant_clock_and_nec_domain_guards(kwargs):
    from p8_determinant import background
    with pytest.raises((TypeError, ValueError)):
        background.reconstruct(**({"Gs": [1], "ys": [1], "cs": [1], "H": 0, "nulls": [0]} | kwargs))
