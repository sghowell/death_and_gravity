"""Modelling scope of Theorem A: exact-identity tests for p4.derivation (sympy + Fractions).

Run:  uv run --with sympy pytest problems/P4/tests/test_derivation.py -q
"""
from fractions import Fraction as Fr

import pytest

sp = pytest.importorskip("sympy")

from p4.derivation import bianchi
from p4.derivation import einstein_euler as ee
from p4.derivation import independent_check as ic
from p4.derivation.qjet import QJet


@pytest.fixture(scope="module")
def D():
    return ee.derive()


# --- deliverable 1: symbolic derivation ---------------------------------------------------------
def test_derivation_sanity(D):
    assert D["u_norm"] == -1
    assert D["U_angular"] == (0, 0)
    assert all(v == 0 for v in D["E_offdiag"])
    assert not D["jets_tr"]["U_t"].has(ee.al_t) and not D["jets_tr"]["U_r"].has(ee.al_t)


def test_textbook_constraint_forms(D):
    """G_tt, G_tr, G_rr of the polar-areal metric in the standard closed forms (exact)."""
    G = D["Einstein"]
    al, aa, r = ee.al, ee.aa, ee.r
    tt = al**2 / (aa**2 * r**2) * (2 * r * ee.aa_r / aa + aa**2 - 1)
    tr = 2 * ee.aa_t / (r * aa)
    rr = (2 * r * ee.al_r / al + 1 - aa**2) / r**2
    assert sp.cancel(ee.to_jet(G[0, 0]) - tt) == 0
    assert sp.cancel(ee.to_jet(G[0, 1]) - tr) == 0
    assert sp.cancel(ee.to_jet(G[1, 1]) - rr) == 0


def test_schwarzschild_is_vacuum(D):
    M = sp.Symbol("M", positive=True)
    sub = {ee.alpha: sp.sqrt(1 - 2 * M / ee.r), ee.a: 1 / sp.sqrt(1 - 2 * M / ee.r)}
    R = D["Ricci"].subs(sub).doit()
    assert all(sp.simplify(R[i, j]) == 0 for i in range(4) for j in range(4))


def test_kha_rows_are_exact_reduction(D):
    out = ee.identities(D)
    assert all(v == 0 for v in out["residuals"].values()), out["residuals"]
    N, a_, W, V, es, ex = ee.N_, ee.a_, ee.W_, ee.V_, ee.es, ee.ex
    detM = W**2 * es**6 / (48 * sp.pi**2 * N**3 * a_**8 * ex**7 * (1 - V**2))
    assert sp.cancel(out["M"].det() - detM) == 0
    assert out["M"][1, 0] == 0


def test_certified_polynomial_system_is_the_derived_system():
    pytest.importorskip("flint")
    assert ee.identities_with_certified_system() == [0, 0, 0, 0]


def test_bianchi_and_angular_equation():
    assert bianchi.contracted_bianchi() == [0, 0, 0, 0]
    # non-vacuity: the divergence of R^{mn} (not G^{mn}) does not vanish
    g, ginv = ee.metric()
    Gam = ee.christoffel(g, ginv)
    R_up = (ginv * ee.ricci(Gam) * ginv).applyfunc(sp.cancel)
    assert any(v != 0 for v in ee.divergence(R_up, Gam))
    c, dt = bianchi.angular_equation_coefficient()
    assert sp.cancel(c - bianchi.EXPECTED_COEFFICIENT) == 0 and dt == 0
    assert bianchi.angular_structure() == 0


# --- deliverable 2: independent exact pointwise check --------------------------------------------
def test_qjet_exact_arithmetic():
    x = QJet.var(2, 2, 0, Fr(3, 2))
    y = QJet.var(2, 2, 1, Fr(-1, 3))
    f = (x * x + y) / (1 + x * y)
    g = f * (1 + x * y) - (x * x + y)
    assert g == 0
    assert (1 / x).d(0).val == -1 / Fr(3, 2) ** 2
    assert f.d(0).d(1).val == f.d(1).d(0).val


def test_independent_pointwise_checks():
    reports = ic.run(n_points=10, seed=20260829)
    assert len(reports) == 10
    assert all(v is True for rep in reports for k, v in rep.items() if k != "M")
    assert ic.cross_check_css() < 1e-9


def _rat(fr):
    return sp.Rational(fr.numerator, fr.denominator)


def test_pointwise_M_and_rows_match_symbolic_exactly():
    """Ties (1) and (2): at a rational point, the symbolic M and rows of (1) evaluate to the
    Fractions of (2) exactly (4 pi absorbed: U_pointwise = 4 pi U_symbolic)."""
    import random
    rng = random.Random(7)
    _, rows_sym = ee.kha_rows_symbolic()
    for _ in range(3):
        p = ic.constrain(ic.random_point(rng))
        k = ic.kha_jets(p)
        a0 = p.a.val
        sub = {ee.a_: _rat(a0), ee.a_s: _rat(k["As"] / (2 * a0)), ee.a_x: _rat(k["Ax"] / (2 * a0)),
               ee.N_: _rat(k["N"]), ee.N_s: _rat(k["Ns"]), ee.N_x: _rat(k["Nx"]),
               ee.W_: _rat(k["W"]), ee.W_s: _rat(k["Ws"]), ee.W_x: _rat(k["Wx"]),
               ee.V_: _rat(k["V"]), ee.V_s: _rat(k["Vs"]), ee.V_x: _rat(k["Vx"]),
               ee.es: _rat(k["es"]), ee.ex: _rat(k["ex"])}
        rw = ic.kha_rows(k)
        for name in ("row1", "row2", "rowM", "row3", "row4"):
            assert rows_sym[name].subs(sub) == _rat(rw[name]), name
        ok, M_pt, _ = ic.euler_vs_rows(p)
        assert ok
        M_sym = (4 * sp.pi * ee.identities()["M"]).subs(sub)
        for i in range(2):
            for j in range(2):
                assert M_sym[i, j] == _rat(M_pt[i][j]), (i, j)
