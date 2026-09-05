import pytest
import sympy as sp
from p8_uv import vacuum as v


def test_vacuum_known_answers():
    assert set(v.derive()["exact_residuals"].values()) == {"0"}
    assert v.forward_coefficient(v.contact()) == 4*v.coupling


def test_direct_component_plane_waves_independent_of_dot_matrix():
    labels = sp.symbols("e0:4")
    # m=3; two incoming and two outgoing particles; sum k_i=0.
    momenta = ((5, 0, 0, 4), (5, 0, 0, -4), (-5, -4, 0, 0), (-5, 4, 0, 0))
    gradients = [sp.I*sum(momenta[j][mu]*labels[j] for j in range(4)) for mu in range(4)]
    kinetic = gradients[0]**2-sum(component**2 for component in gradients[1:])
    expanded = sp.Poly(sp.expand(v.coupling*kinetic**2), *labels).coeff_monomial(sp.prod(labels))
    assert expanded == v.contact().subs({v.mass: 3, v.s: 100, v.transfer: -32, v.w: -32})


def test_actual_heavy_threshold_and_mass_correction():
    b2 = v.forward_coefficient(v.heavy_exchange())
    assert b2.subs({v.mass: 1, v.heavy_mass: 3, v.mediator: 1}) == sp.Rational(1, 14)
    leading = v.mediator**2/(2*v.heavy_mass**2)
    assert sp.cancel(b2-leading-v.mediator**2*v.mass**2/
                     (v.heavy_mass**2*(v.heavy_mass**2-2*v.mass**2))) == 0


def test_bad_dot_matrix_rejected():
    with pytest.raises(ValueError, match="four-leg"):
        v.labelled_contact(sp.eye(3))
    wrong = sp.eye(4)
    wrong[0, 1] = 1
    with pytest.raises(ValueError, match="four-leg"):
        v.labelled_contact(wrong)
