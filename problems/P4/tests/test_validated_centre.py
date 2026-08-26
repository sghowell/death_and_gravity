"""S2 tests for the validated regular-centre expansion (Theorem A, item A2)."""
import numpy as np
import pytest
from flint import arb
from p4 import shoot
from p4.validated import centre
from p4.validated.arbseries import precision

NH, WH = "1.2365999612", "5.82098013"          # S1 limits (N^_inf, w^_inf); V^_inf = -1/(2 N^_inf)
V0_EC = 0.112439401388092


@pytest.fixture(scope="module")
def ex():
    return centre.centre_expansion(WH, nhat=NH, K=30)


def test_fixed_point_and_spectrum(ex):
    n0, w0, v0 = ex.balls[0]
    # N^ V^ = -1/2 gives V^_inf = -0.4043344782; S1's quoted -0.40433391 carries its e^{-3x} drift (5.7e-7)
    assert abs(float(v0) - (-0.40433391)) < 1e-6
    assert (n0 * v0 + arb(1) / 2).contains(arb(0))
    D = np.array([[float(ex.D[i, j]) for j in range(3)] for i in range(3)])
    E = np.array([[float(ex.E[i, j]) for j in range(3)] for i in range(3)])
    eig = np.sort(np.linalg.eigvals(-np.linalg.solve(D, E)).real)
    assert np.allclose(eig, [-3, 0, 0], atol=1e-10)                    # exponents {-3, 0, 0}: no resonance
    fl = ex.floats()
    assert np.abs(fl[:, 1::2]).max() < 1e-60                            # only even powers of e^x


def test_reproduces_S1_profile(ex):
    """(iv) with S1's parameters the series reproduces S1's integrated profile near the centre."""
    sh = shoot.shoot(V0_EC, 1, x_end=-10.0, keep_sol=True)
    for x, tol in ((-4.0, 1e-7), (-6.0, 1e-7)):
        ref = sh.sol.sol(x)                                            # (N^, w^, V^) from S1
        with precision(256):
            vals = [float(v) for v in ex.eval(x, with_tail=False)]
        assert np.allclose(vals, ref, atol=tol), (x, vals, ref)
    assert abs(float(ex.Ahat_series()(arb("0.0025"))) - 2 * float(ex.balls[0][1]) / 3) < 1e-3


def test_translation_scaling_identity(ex):
    with precision(256):
        nh, wh = arb(NH), arb(WH)
        exn = centre.centre_expansion(nh**2 * wh, nhat=1, K=30)
        resc = exn.rescale(nh)
        for k in range(31):
            for i in range(3):
                assert (resc[k][i] - ex.balls[k][i]).contains(arb(0))


def test_tail_certificate_valid_to_x_minus_3(ex):
    cert = ex.certify()
    assert cert.ok and float(cert.nu) >= np.exp(-3.0)
    with precision(256):
        big = centre.centre_expansion(WH, nhat=NH, K=45)
        nu, eps = cert.nu, cert.eps
        tail = arb(0)
        for k in range(31, 46):
            uk = max(abs(c) for c in big.balls[k])
            assert uk <= eps * nu ** (-k)
            tail += uk * nu**k
        assert tail <= eps
        assert float(cert.tail_bound(arb(-6.0).exp().abs_upper())) < 1e-20


def test_parametrised_family_taylor_model():
    with precision(256):
        mu = arb(NH) ** 2 * arb(WH)
    exi = centre.centre_expansion(mu, nhat=1, K=20, width=1e-8, m=3)
    rad = exi.radii()
    fl = np.abs(exi.floats())
    assert rad[:, 0].max() < 2e-8
    assert (rad[:, 1:] <= 1e-8 * 20 * np.arange(1, 21) * (fl[:, 1:].max(axis=0) + 1)).all()  # ~ k|Y_k| w
    assert max(max(r) for r in exi.rem) < 1e-10
    with precision(256):
        exp_pt = centre.centre_expansion(mu, nhat=1, K=20)
        for k in range(21):
            for i in range(3):
                assert exi.balls[k][i].contains(exp_pt.balls[k][i])
