import sympy as sp
from p8_scattering import amplitude, frequency, verify


def test_three_exchange_partitions_keep_scalar_and_both_tensor_polarizations():
    out = verify.frozen_example("ssss", 10, 0, "gamma")
    assert len(out['channels']) == 3
    assert all([e['kind'] for e in c['entries']] == ['s', 't', 't'] for c in out['channels'])
    assert any(e['contribution'] != '0' for c in out['channels'] for e in c['entries'] if e['kind'] == 't')
    assert out['external_shell_residuals'] == ['0']*4
    assert out['contact'] != out['total']


def test_known_answer_free_redefinition_from_labelled_lagrangian_jets():
    E, g = sp.symbols('E g', real=True)
    waves = ((0, 3, 4), (0, 3, -4), (-4, -3, 0), (4, -3, 0))
    energies = (E, E, -E, -E)
    def vertex(selected, momenta, energy):
        c = frequency.Context(momenta, (E, g))
        field = sum(c.leg(i) for i in range(len(momenta)))
        velocity = sum(-sp.I*w*c.leg(i) for i, w in enumerate(energy))
        kinetic = velocity**2-sum(field.derivative(i)**2 for i in range(3))
        expression = g*field*kinetic-g*field**3/2 if selected == 3 else g*g*field**2*kinetic/2-g*g*field**4/8
        return sp.expand(expression.coefficient(c.full))
    contact = vertex(4, waves, energies).subs(E**2, 26)
    total = contact
    for left, right in frequency.PARTITIONS:
        k = tuple(-sum(waves[i][j] for i in left) for j in range(3))
        w = -sum(energies[i] for i in left)
        v1 = vertex(3, tuple(waves[i] for i in left)+(k,), tuple(energies[i] for i in left)+(w,))
        v2 = vertex(3, tuple(waves[i] for i in right)+(tuple(-v for v in k),), tuple(energies[i] for i in right)+(-w,))
        inverse = w*w-sum(v*v for v in k)-1
        total -= sp.cancel(v1*v2/inverse).subs(E**2, 26)
    assert contact == g*g
    assert sp.factor(total) == 0


def test_real_action_frequency_reversal():
    legs = verify.fixture("stst", 10)[:2]
    # Use a cubic closed triangle with real rational frequencies.
    left = frequency.Leg(legs[0].wave, 's', 3)
    middle = frequency.Leg(legs[1].wave, 't', 7, legs[1].polarization)
    wave = tuple(-left.wave[i]-middle.wave[i] for i in range(3))
    inputs = (left, middle, frequency.Leg(wave, 's', -10))
    value = frequency.vertex(inputs, sp.Rational(3, 5), 'unitary')['kernel']
    reversed_value = frequency.vertex(tuple(frequency.reverse(l) for l in inputs), sp.Rational(3, 5), 'unitary')['kernel']
    assert sp.simplify(reversed_value-sp.conjugate(value)) == 0


def test_mixed_and_tensor_external_shells():
    for kinds in ('stst', 'tttt'):
        out = verify.frozen_example(kinds, 10, 0, 'gamma')
        assert out['external_shell_residuals'] == ['0']*4
        assert out['total'] != '0'


def test_frozen_example_is_not_silently_called_a_physical_s_matrix():
    out = amplitude.tree(verify.fixture('ssss', 10), 0)
    assert 'not a time-dependent physical amplitude' in out['scope']
