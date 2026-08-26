"""Stage C tests: kappa_1, the gauge eigenvalue, sonic-point perturbation data."""
import numpy as np
import pytest

from p4 import css, perturb, shoot, spectrum, taylor

V0_EC = 0.112439401388092


@pytest.fixture(scope="module")
def prob():
    return perturb.Problem(V0_EC, 1, K=36, x_end=-9.0)


def test_perturbation_series_satisfies_linearized_constraint(prob):
    for kap in (2.81, 0.3557, 1.5 + 2.0j):
        pc, sv, cons = taylor.perturbation_series(prob.bg, kap, K=24, return_info=True)
        assert pc[0, 0] == 1.0 and pc[1, 0] == 0.0                # normalisation, gauge
        assert np.abs(cons[:22]).max() < 1e-9                     # kappa A_p = dC . y_p
        assert np.nanmin(sv[1:]) > 1e-3                           # no resonance


def test_resonances_outside_the_box(prob):
    """Second Frobenius exponent sigma(kappa) = n only at kappa = -0.099 - 1.099 n."""
    for n in (1, 2, 3):
        kr = taylor._det_root(prob.bg, n)
        assert kr.real < -1.0 and abs(kr.imag) < 1e-8
        assert abs(kr.real - (-0.0990 - 1.0990 * n)) < 2e-3


def test_kappa_1(prob):
    k1, absE, it = spectrum.refine_zero(prob, 2.81, 2.812, real=True)
    assert abs(k1 - 2.81055255) < 1e-6                            # KHA99 Table 2
    assert absE < 1e-10


def test_gauge_eigenvalue_found_and_identified(prob):
    kg, absE, it = spectrum.refine_zero(prob, 0.355, 0.357, real=True)
    assert abs(kg - shoot.gauge_kappa(V0_EC)) < 1e-8              # = -dN_ss/dx(0)
    # sonic-point data of the matching solution equals the pure-gauge generator
    y0, br = css.sonic_branches(V0_EC)
    y1 = br[1]
    gauge = np.array([y1[0], y1[1] + kg * y0[1], y1[2], y1[3]]) / y1[0]
    pc = taylor.perturbation_series(prob.bg, kg, K=6)
    assert np.allclose(pc[:, 0].real, gauge, atol=1e-8)


def test_centre_has_one_irregular_direction(prob):
    for kap in (0.5, 2.81, 3.0 + 4.0j):
        ev, _ = perturb.centre_exponents(prob, kap)
        ev = np.sort(ev.real)
        assert abs(ev[0] + 3.0) < 1e-4                            # mass mode, A_p ~ e^{-x}
        assert np.all(np.abs(ev[1:]) < 1e-4)                      # two regular directions


def test_winding_number_small_box(prob):
    """Exactly two zeros (gauge + kappa_1) in [0, 4] x [-1, 1]."""
    w, ks, Es = spectrum.winding_number(prob, (0.0, 4.0, -1.0, 1.0), n_per_side=24)
    assert w == 2
