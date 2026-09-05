import pytest
import sympy as sp
from p8_physical import lagrangian
from p8_physical.vertices import Leg, hamiltonian_kernel, tensor_basis
from p8_s5 import compact

K = (10, 20, 30)
OPPOSITE = tuple(-v for v in K)
Q = sum(v*v for v in K)


@pytest.mark.parametrize("kinds", [("s", "s"), ("s", "p"), ("p", "p")])
def test_full_scalar_quadratic_matches_pinned_hamiltonian(kinds):
    x, J = compact.x, 2*compact.P
    theta, lam = 2*x, 1-2*(1-x*x)**3
    expected = {("s", "s"): 2*(lam**2*Q**2/J-Q),
                ("s", "p"): -theta*lam*Q/J, ("p", "p"): theta**2/(2*J)}
    actual = hamiltonian_kernel((Leg(K, kinds[0]), Leg(OPPOSITE, kinds[1])))["kernel"]
    assert sp.cancel(actual-expected[kinds]) == 0


@pytest.mark.parametrize("kinds", [("s", "s"), ("s", "p"), ("p", "p")])
def test_gamma_canonical_swap_includes_time_dependent_generator(kinds):
    x, J = compact.x, 2*compact.P
    theta, lam = 2*x, 1-2*(1-x*x)**3
    expected = {("s", "s"): 2*theta**2*Q**2/J,
                ("s", "p"): theta*lam*Q/J-2*x,
                ("p", "p"): (Q*lam**2-J)/(2*J*Q)}
    actual = hamiltonian_kernel((Leg(K, kinds[0]), Leg(OPPOSITE, kinds[1])), chart="gamma")["kernel"]
    assert sp.cancel(actual-expected[kinds]) == 0
    if kinds == ("s", "p"):
        assert sp.cancel(actual-theta*lam*Q/J) != 0  # missing generator is detected


@pytest.mark.parametrize("kinds", [("t", "t"), ("pi", "pi"), ("t", "pi"), ("s", "t")])
def test_tensor_kinetic_gradient_and_scalar_tensor_decoupling(kinds):
    E = tensor_basis(K)[0]
    norm = sp.trace(sp.Matrix(E)**2)
    legs = tuple(Leg(k, kind, E if kind in ("t", "pi") else None)
                 for k, kind in zip((K, OPPOSITE), kinds))
    expected = {("t", "t"): Q*norm/4, ("pi", "pi"): 4/norm, ("t", "pi"): 0, ("s", "t"): 0}
    assert sp.cancel(hamiltonian_kernel(legs)["kernel"]-expected[kinds]) == 0


def test_background_canonical_tadpoles_vanish():
    for kind in ("s", "p"):
        assert hamiltonian_kernel((Leg((0, 0, 0), kind),))["kernel"] == 0


@pytest.mark.parametrize("chart,point", [("unitary", sp.Rational(3, 5)), ("gamma", 0)])
@pytest.mark.parametrize("kind", ["s_dot", "t_dot"])
def test_exact_unit_kinetic_normalization(chart, point, kind):
    E = tensor_basis(K)[1] if kind == "t_dot" else None
    legs = (Leg(K, kind, E), Leg(OPPOSITE, kind, E))
    assert lagrangian.kernel(legs, point, chart)["kernel"] == 1
