"""Physical and source dictionaries without unsupported observational input."""

import sympy as sp
from p8a_maxwell_cosmology import dictionary


def test_every_proper_Planck_and_perfect_fluid_identity():
    assert all(sp.simplify(value) == 0 for value in dictionary.identities().values())


def test_actual_observer_scale_only_equals_reference_if_anchored():
    free = dictionary.proper_scales(10)
    fixed = dictionary.proper_scales(10, anchored=True)
    assert free["reference_Hstar"] == sp.Rational(1, 5)
    assert free["actual_observed_H0_range"] == [sp.Rational(199, 1000), sp.Rational(201, 1000)]
    assert fixed["actual_observed_H0_range"] == [sp.Rational(1, 5)]*2


def test_negative_extra_EED_and_positive_Lambda_cannot_be_omitted():
    bad = dictionary.source_budget(kappa=2, hbar=0, tau=1, beta_m=0,
                                    cosmological_constant=1, other_eed_lower=-3)
    good = dictionary.source_budget(kappa=2, hbar=0, tau=1, beta_m=0,
                                     cosmological_constant=-1, other_eed_lower=3)
    assert bad["sigma"] == 7 and bad["sufficient_source_budget"] is False
    assert good["sigma"] == 0 and good["sufficient_source_budget"] is True


def test_finite_beta_is_retained_in_the_physical_joint_budget():
    data = dictionary.source_budget(kappa=1, hbar=1, tau=10000, beta_m=-71,
                                     cosmological_constant=0, other_eed_lower=0)
    assert data["rational_weighted_delta_upper"] == sp.Rational(1, 10**8)
    assert data["sufficient_quantum_budget"] is True
    assert data["beta_M"] == -71
