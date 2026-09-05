import pytest
import sympy as sp
from p8_physical import lagrangian, momentum
from p8_physical.vertices import Leg, chart_functions, hamiltonian_kernel, tensor_basis
from p8_s5 import compact
from p8_s5 import nonlinear_d as m

QUAD = ((10, 20, 0), (0, 10, 30), (20, -10, 10), (-30, -20, -40))
TRI = ((10, 20, 0), (0, 10, 30), (-10, -30, -30))


def test_quartic_legendre_identity_and_partition_factor():
    assert lagrangian.stationary_checks()["quartic_Legendre"] == 0
    # H3=g*p*q^2, H2=p^2/2 gives L4=g^2*q^4/2: labelled coefficient 12*g^2.
    g = sp.Symbol("g")
    assert sum((2*g)**2 for _ in range(3)) == sp.factorial(4)*g**2/2


def test_mixed_tensor_legendre_channels_are_not_dropped():
    kinds = ("s", "s", "s_dot", "s_dot")
    legs = tuple(Leg(k, kind) for k, kind in zip(QUAD, kinds))
    out = lagrangian.kernel(legs, 0, "gamma", False)
    assert any(v["rational_tensor"] != 0 for v in out["Legendre_corrections"])
    assert sp.cancel(out["kernel"]-out["rational_minus_H"]) != 0
    assert len(out["Legendre_corrections"]) == 3


def test_permutation_of_external_legs_preserves_quartic_kernel():
    legs = tuple(Leg(k, "t" if i == 3 else "s", tensor_basis(k)[0] if i == 3 else None)
                 for i, k in enumerate(QUAD))
    original = hamiltonian_kernel(legs, sp.Rational(3, 5))["kernel"]
    assert original != 0
    permuted = (legs[2], legs[0], legs[3], legs[1])
    assert hamiltonian_kernel(permuted, sp.Rational(3, 5))["kernel"] == original


def test_rotation_of_all_momenta_and_polarizations_preserves_mixed_kernel():
    rotation = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
    legs = (Leg(TRI[0], "s"), Leg(TRI[1], "pi", tensor_basis(TRI[1])[0]),
            Leg(TRI[2], "pi", tensor_basis(TRI[2])[1]))
    transformed = tuple(Leg(tuple(rotation*sp.Matrix(l.wave)), l.kind,
                            tuple(map(tuple, (rotation*sp.Matrix(l.polarization)*rotation.T).tolist()))
                            if l.polarization else None) for l in legs)
    assert sp.cancel(hamiltonian_kernel(legs)["kernel"]-hamiltonian_kernel(transformed)["kernel"]) == 0


def test_fixed_local_unit_derivative_includes_normalization_weight():
    u, d = m.u, m.d
    q = sp.Symbol("q", positive=True)
    k2 = sp.Symbol("k2", positive=True)
    assert sp.simplify(sp.sqrt(d)*sp.diff(u/sp.sqrt(d), u)-1/d) == 0
    assert sp.simplify(sp.sqrt(d)*sp.diff(k2/d, u)+2*u/sp.sqrt(d)*(k2/d)) == 0
    f = m.functions()
    K_old = (k2/d**2)*f["J"]/((k2/d**2)*f["Lambda"]**2-f["J"])
    actual = sp.cancel(sp.diff(K_old, u)/(2*K_old))
    point = sp.Rational(3, 4)
    expected = chart_functions("gamma", q)["logZdot"].subs({compact.x: sp.Rational(3, 5), q: 100})
    evaluated = actual.subs({u: point, k2: sp.Rational(25, 16)*100})*sp.Rational(5, 4)
    assert sp.cancel(evaluated-expected) == 0


def test_time_dependent_normalization_is_more_than_dividing_external_legs():
    legs = tuple(Leg(k, "s") for k in TRI)
    normal = lagrangian.kernel(legs, sp.Rational(3, 5), "unitary", True)
    raw = lagrangian.kernel(legs, sp.Rational(3, 5), "unitary", False)
    assert sp.cancel(normal["rational_kernel"]-raw["rational_kernel"]) != 0


def test_exceptional_wavevector_and_bad_velocity_chart_fail_closed():
    legs = tuple(Leg(k, "s") for k in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)))
    with pytest.raises(momentum.ZeroMomentumConstraint):
        hamiltonian_kernel(legs, 0)
    pair = (Leg((10, 0, 0), "s_dot"), Leg((-10, 0, 0), "s_dot"))
    with pytest.raises(ValueError, match="Singular"):
        lagrangian.kernel(pair, 0, "unitary")
    soft_internal = tuple(Leg(k, "s_dot") for k in ((10, 0, 0), (-9, 0, 0), (0, 10, 0), (-1, -10, 0)))
    with pytest.raises(ValueError, match="Internal momentum"):
        lagrangian.kernel(soft_internal, 0, "gamma")
