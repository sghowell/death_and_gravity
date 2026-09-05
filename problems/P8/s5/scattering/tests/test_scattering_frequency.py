from itertools import product

import pytest
import sympy as sp
from p8_physical import lagrangian
from p8_physical.vertices import Leg as OldLeg
from p8_physical.vertices import tensor_basis
from p8_scattering import amplitude, frequency

TRI = ((10, 20, 0), (0, 10, 30), (-10, -30, -30))
QUAD = ((10, 20, 0), (0, 10, 30), (20, -10, 10), (-30, -20, -40))


@pytest.mark.parametrize("kinds,chart,x", [(('s', 's', 's'), 'gamma', 0),
                                          (('s', 't', 't'), 'unitary', sp.Rational(3, 5)),
                                          (('s', 's', 's', 's'), 'gamma', 0),
                                          (('s', 't', 's', 't'), 'unitary', sp.Rational(3, 5)),
                                          (('t', 't', 't', 't'), 'gamma', 0)])
def test_direct_frequency_seeds_equal_independent_velocity_multilinearity(kinds, chart, x):
    waves = TRI if len(kinds) == 3 else QUAD
    legs = tuple(frequency.Leg(k, kind, i+1, tensor_basis(k)[0] if kind == 't' else None)
                 for i, (k, kind) in enumerate(zip(waves, kinds)))
    actual = frequency.vertex(legs, x, chart)["kernel"]
    expected = 0
    for flags in product((False, True), repeat=len(legs)):
        inputs = tuple(OldLeg(l.wave, l.kind+'_dot' if flag else l.kind, l.polarization)
                       for l, flag in zip(legs, flags))
        factor = sp.prod(-sp.I*l.frequency-3*x for l, flag in zip(legs, flags) if flag)
        expected += factor*lagrangian.kernel(inputs, x, chart)["kernel"]
    assert sp.simplify(actual-expected) == 0


def test_free_redefinition_requires_both_contact_and_exchange():
    out = amplitude.free_redefinition_check()
    assert out['residual'] == 0
    assert out['contact'] != 0
    assert sum(out['exchange']) != 0


def test_frequency_total_is_enforced():
    legs = tuple(frequency.Leg(k, 's', 1) for k in QUAD)
    with pytest.raises(ValueError, match="frequency zero"):
        amplitude.tree(legs, 0)
