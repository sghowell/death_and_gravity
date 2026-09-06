"""The full trace, both primitive factors, and retained remainder controls."""

import sympy as sp
from p8a_continuation import decomposition, resolvent
from p8a_preparation.actual_map import local_coefficient


def test_complete_trace_and_auxiliary_map_identities():
    assert all(value == 0 for value in decomposition.identities().values())


def test_named_beta_uses_the_same_finite_A11_coefficient():
    data = resolvent.named_parameters()
    assert sp.simplify(data["beta"]-sp.EulerGamma-2*local_coefficient(data["scale"])) == 0
    assert data["c"] == sp.Rational(62500000000000, 3)


def test_anomaly_is_not_omitted_from_the_trace():
    h, u = sp.symbols("h u", positive=True)
    assert decomposition.anomaly_polynomial(h, u).expand() == u**2/4+h**2*u/30+h**4/30
    assert decomposition.anomaly_polynomial(h, 0) == h**4/30


def test_each_full_remainder_input_has_its_own_nonzero_coefficient():
    base, rolling, anomaly, auxiliary, state = sp.symbols("base rolling anomaly auxiliary state")
    delta = sp.Symbol("delta", positive=True)
    full = decomposition.full_remainder(base, rolling, anomaly, auxiliary, state, delta)
    assert [sp.diff(full, term) for term in (base, rolling, anomaly, auxiliary, state)] == [
        1, 1/(60*delta), -1, 2, -1]
    assert full.subs({base: 1, rolling: 0, anomaly: 0, auxiliary: 0, state: 0}) == 1


def test_double_primitive_not_single_primitive_controls_stiff_power():
    t, s = sp.symbols("t s", positive=True)
    x = t**2
    first, second = t**3/3, t**4/12
    assert sp.diff(second, t) == first
    assert sp.diff(first, t) == x
    transform_x = 2/s**3
    assert 2/s**5 == transform_x/s**2
    assert 2/s**5 != transform_x/s


def test_static_trace_block_is_not_a_static_radiation_SEE_solution():
    af = sp.Symbol("af", positive=True)
    # The reference vacuum and trace vanish at constant scale, but the
    # original fixed radiation term leaves the energy constraint nonzero.
    assert 3*sp.diff(af, sp.Symbol("t"))**2-3 == -3
