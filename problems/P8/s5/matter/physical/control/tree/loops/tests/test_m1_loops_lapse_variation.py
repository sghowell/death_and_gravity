"""Independent lapse-first variation, not a conservation-derived pressure oracle."""

from functools import cache

import sympy as sp


@cache
def variation():
    t = sp.Symbol("t", real=True)
    a, lapse = sp.Function("a")(t), sp.Function("N")(t)
    ricci = -6*(sp.diff(a, t, 2)/(a*lapse**2)
                + sp.diff(a, t)**2/(a**2*lapse**2)
                - sp.diff(a, t)*sp.diff(lapse, t)/(a*lapse**3))
    action = lapse*a**3*ricci**2

    def euler_lagrange(field, order, lagrangian):
        return sum((-1)**j*sp.diff(sp.diff(lagrangian, sp.diff(field, t, j)), t, j)
                   for j in range(order+1))

    gauge = {lapse: 1, **{sp.diff(lapse, t, j): 0 for j in range(1, 5)}}
    density = sp.factor((-lapse**2*euler_lagrange(lapse, 1, action)/(2*a**3)).subs(gauge))
    pressure = sp.factor((euler_lagrange(a, 2, action)/(6*lapse*a**2)).subs(gauge))
    wrong_action = action.subs(sp.diff(lapse, t), 0)
    wrong_density = sp.factor(
        (-lapse**2*euler_lagrange(lapse, 1, wrong_action)/(2*a**3)).subs(gauge))
    return t, a, density, pressure, wrong_density


def specialize(expression, a, t, value):
    jets = {sp.diff(a, t, j): sp.diff(value, t, j) for j in range(5)}
    return sp.factor(expression.subs(jets, simultaneous=True))


def test_lapse_and_scale_variations_match_both_covariant_components():
    t, a, density, pressure, _ = variation()
    hubble = sp.diff(a, t)/a
    ricci = -6*(sp.diff(hubble, t)+2*hubble**2)
    covariant_density = (-6*ricci*(sp.diff(hubble, t)+hubble**2)-ricci**2/2
                         + 6*hubble*sp.diff(ricci, t))
    covariant_pressure = (2*ricci*(sp.diff(hubble, t)+3*hubble**2)+ricci**2/2
                          - 2*sp.diff(ricci, t, 2)-4*hubble*sp.diff(ricci, t))
    assert sp.factor(density-covariant_density) == 0
    assert sp.factor(pressure-covariant_pressure) == 0
    assert sp.factor(sp.diff(density, t)+3*hubble*(density+pressure)) == 0


def test_cd_components_from_independent_lapse_first_variation():
    t, a, density, pressure, _ = variation()
    tau = sp.Symbol("tau", positive=True)
    u = t/tau
    cd = (1+u**2)**2
    assert sp.factor(specialize(density, a, t, cd)
                     - 288*(21*u**4-14*u**2+1)/(tau**4*(1+u**2)**4)) == 0
    assert sp.factor(specialize(pressure, a, t, cd)
                     + 576*(7*u**4-1)/(tau**4*(1+u**2)**4)) == 0


def test_radiation_variation_vanishes_not_just_its_ricci_scalar():
    t, a, density, pressure, _ = variation()
    assert specialize(density, a, t, sp.sqrt(t)) == 0
    assert specialize(pressure, a, t, sp.sqrt(t)) == 0


def test_omitting_lapse_derivative_before_variation_fails_at_bounce():
    t, a, density, _, wrong_density = variation()
    mismatch = specialize(wrong_density-density, a, t, (1+t**2)**2)
    assert mismatch.subs(t, 0) == 576
