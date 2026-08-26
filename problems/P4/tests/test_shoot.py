"""Stage B tests: the Evans-Coleman root, its regular centre, zero count, gauge value."""
import numpy as np
import pytest

from p4 import css, shoot

V0_EC = 0.112439401388092


@pytest.fixture(scope="module")
def ec_root():
    V0, sh = shoot.refine_root(0.105, 0.120, 1, x_end=-10.0)
    return V0, sh


def test_ec_root_value(ec_root):
    V0, sh = ec_root
    assert abs(V0 - V0_EC) < 1e-10
    assert sh.reason == "end"


def test_ec_regular_centre(ec_root):
    V0, sh = ec_root
    # V0 is determined to ~1e-13 (noise floor of the mismatch); the residual mass mode
    # grows like e^{-3x}, so test at x = -6 where the profile is regular to ~1e-5.
    x = -6.0
    yh = css.scaled3_to_full(x, sh.sol.sol(x))
    cc = css.centre_check(x, yh)
    assert abs(cc["Ah_minus_2Wh3"]) < 1e-3            # A_{-inf} = 2 w_{-inf}/3   (KHA99 222)
    assert abs(cc["NV_plus_half"]) < 1e-4             # N_{-inf} V_{-inf} = -1/2
    assert abs(cc["mass"]) < 1e-6
    assert yh[0] > 0 and yh[2] > 0                    # A >= 1, omega > 0


def test_ec_single_zero_of_V(ec_root):
    V0, sh = ec_root
    assert sh.n_zeros_V == 1
    n_in, n_out, inner, outer = shoot.count_zeros_V_full(V0, 1, x_end=-8.0, x_max=20.0)
    assert n_in == 1 and n_out == 0
    assert outer.status == 0


def test_ec_t0_limits(ec_root):
    V0, _ = ec_root
    out = shoot.continue_outward(V0, 1, x_max=25.0)
    N, W, V = out.y[:, -1]
    Am1 = css.Am1_constraint(N, W, V)
    a, m_over_r = np.sqrt(1 + Am1), Am1 / (1 + Am1) / 2
    assert abs(a - 1.0653) < 2e-3            # EC94: 1.07
    assert abs(m_over_r - 0.05945) < 5e-4    # EC94: 0.0596


def test_gauge_kappa(ec_root):
    V0, _ = ec_root
    kb = shoot.gauge_kappa(V0)
    assert abs(kb - 0.3556992) < 1e-6        # KHA99 sec. V.7.2 quotes 0.35699 (see report)


def test_other_branch_has_no_root():
    shots = shoot.scan(np.linspace(0.02, 0.55, 12), 0, x_end=-6.0)
    assert shoot.brackets_from_scan(shots) == []
