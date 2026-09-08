"""Continuous proof arithmetic, independent jets, preparation and unit checks."""
import pytest
import sympy as sp
from p8_affine_aligned import modes
from p8_affine_retuned.bounds import units
from p8_aligned_response import bounds as b
from p8_aligned_response import derivatives as d


def test_exact_derivative_energy_and_quasistatic_identities():
    assert all(value == 0 for value in d.checks().values())
    assert all(value == 0 for value in b.checks().values())
    assert all(value is True for value in d.proof_checks().values())
    assert all(value is True for value in b.proof_checks().values())


@pytest.mark.parametrize("point", (-10, -1, -sp.Rational(1, 2), 0, sp.Rational(1, 4), 1, 10))
def test_actual_frequency_jets_at_fixed_comoving_momentum(point):
    W = modes.canonical()["longitudinal"]["frequency_squared"].subs(
        {modes.ZETA: b.ZETA_MAX, modes.KCOM2: 1})
    assert W.subs(modes.u, point) >= b.D_MIN
    assert abs(sp.diff(W, modes.u).subs(modes.u, point)) <= 78
    assert abs(sp.diff(W, modes.u, 2).subs(modes.u, point)) <= 712


@pytest.mark.parametrize("point", (-sp.Rational(1, 2), -sp.Rational(1, 4), 0, sp.Rational(1, 4), sp.Rational(1, 2)))
def test_prepared_cubic_source_satisfies_canonical_jet_envelope(point):
    data = modes.canonical()
    weight = data["longitudinal"]["weight"].subs({modes.ZETA: b.ZETA_MAX, modes.KCOM2: 1})
    source = (modes.u+sp.Rational(1, 2))**3
    g = sp.sqrt(weight)
    force = sp.diff(weight*source, modes.u)/g
    result = b.bound_values(b.ZETA_MAX, -sp.Rational(1, 2), sp.Rational(1, 2), 1, 6)
    A = result["canonical_force_jet_A"]
    for order in (0, 1, 2):
        actual = sp.factor(sp.diff(force, modes.u, order).subs(modes.u, point))
        assert sp.factor(A**2-actual**2) >= 0
    assert force.subs(modes.u, -sp.Rational(1, 2)) == 0
    assert sp.diff(force, modes.u).subs(modes.u, -sp.Rational(1, 2)) == 0


@pytest.mark.parametrize("point", (-sp.Rational(1, 2), -sp.Rational(1, 4), 0, sp.Rational(1, 4), sp.Rational(1, 2)))
def test_manufactured_exact_retarded_solution_and_error(point):
    # v=(u+1/2)^5 has v=v'=0 and J=J'=0 at the left endpoint.
    # The jet envelope follows directly from |W|<=m²+16, |W'|<=78,
    # |W''|<=712 and the polynomial derivative bounds 1,5,20,60,120.
    v = (modes.u+sp.Rational(1, 2))**5
    W = modes.canonical()["longitudinal"]["frequency_squared"].subs(
        {modes.ZETA: b.ZETA_MAX, modes.KCOM2: 1})
    force = sp.diff(v, modes.u, 2)+W*v
    A = 20/b.ZETA_MAX+1932
    for order in (0, 1, 2):
        assert abs(sp.factor(sp.diff(force, modes.u, order).subs(modes.u, point))) <= A
    assert abs(sp.factor((v-b.ZETA_MAX*force).subs(modes.u, point))) <= b.ZETA_MAX*A/80
    for item in (v, sp.diff(v, modes.u), force, sp.diff(force, modes.u)):
        assert item.subs(modes.u, -sp.Rational(1, 2)) == 0


def test_unprepared_or_nonzero_heavy_data_cannot_be_included():
    # Constant forcing and zero heavy data have v(0)=0 but zeta*J=zeta:
    # the proposed 1/80 error fails at the initial endpoint without J(0)=0.
    assert b.ZETA_MAX > b.ZETA_MAX/80
    # With zero forcing, a unit homogeneous initial displacement has
    # unit error while every source-proportional bound is zero.
    assert 1 > b.bound_values(b.ZETA_MAX, 0, sp.Rational(1, 2), 1, 0)["canonical_error_against_zeta_J"]


def test_homogeneous_scalar_source_is_temporal_and_algebraic():
    result = b.bound_values(b.ZETA_MAX, -sp.Rational(1, 2), sp.Rational(1, 2), 0, 1)
    assert result["branch"] == "homogeneous_algebraic_temporal_source"
    assert result["canonical_force_jet_A"] == result["canonical_error_against_zeta_J"] == 0
    assert result["temporal_readout_error_against_S"] == 0


def test_nonunit_physical_readout_errors():
    conversion = units(3, 2, sp.Rational(1, 10000))
    result = b.bound_values(conversion["normalized_zeta"], -sp.Rational(1, 2), sp.Rational(1, 2), 1, sp.Rational(1, 100))
    assert result["zeta"] == sp.Rational(1, 120000)
    assert result["spatial_readout_error_against_zeta_sqrt_q_times_Sprime_plus_2rhoS"]/2 == sp.Rational(1, 4000000)
    assert result["temporal_readout_error_against_S"]/2 == sp.Rational(1, 48000)
    assert result["full_nonlinear_or_quantum_remainder_claim"] is False


@pytest.mark.parametrize("values", ((0, 0, sp.Rational(1, 2), 1, 1),
                                   (sp.Rational(1, 2000), 0, sp.Rational(1, 2), 1, 1),
                                   (b.ZETA_MAX, -1, 0, 1, 1), (b.ZETA_MAX, 0, 0, 1, 1),
                                   (b.ZETA_MAX, 0, sp.Rational(1, 2), -1, 1),
                                   (b.ZETA_MAX, 0, sp.Rational(1, 2), 2, 1),
                                   (b.ZETA_MAX, 0, sp.Rational(1, 2), 1, -1),
                                   (0.00005, 0, sp.Rational(1, 2), 1, 1),
                                   (True, 0, sp.Rational(1, 2), 1, 1)))
def test_outside_prepared_domain_rejected(values):
    with pytest.raises((TypeError, ValueError)):
        b.bound_values(*values)
