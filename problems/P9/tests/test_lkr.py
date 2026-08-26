"""LKRModel2 (shared row builder) must agree with LKRModel; Verifier3 must reproduce solver dual bounds."""
import numpy as np
import pytest

from p9.bound import sn_subset
from p9.data import load_desi, load_pantheon
from p9.geometry import lcdm_u_nodes
from p9.lkr import LKRModel, initial_brackets3
from p9.lkr2 import LKRModel2
from p9.model import ClassSpec, Frozen
from p9.verify3 import Verifier3


@pytest.fixture(scope="module")
def small():
    bao = load_desi(); sn = load_pantheon(); sub = sn_subset(sn, 40)
    spec = ClassSpec(L=1.5, grid_kind="geometric")
    fr = Frozen(bao, sub, spec)
    u = lcdm_u_nodes(spec.x, 0.30, 101.0)
    T = fr.chi2(u, fr.best_Mp(u)) + 4.0
    br = initial_brackets3(fr)
    return fr, br, T


def test_lkr2_matches_lkr(small):
    fr, br, T = small
    m1 = LKRModel(fr, br, T); m2 = LKRModel2(fr, br, T)
    assert m1.nvar == m2.nvar
    a1 = m1.min_lambda0(); a2 = m2.min_lambda0()
    assert abs(a1 - a2) < 1e-6
    c1, _ = m1.min_chi2(); c2, _ = m2.min_chi2()
    assert abs(c1 - c2) < 1e-4


@pytest.mark.parametrize("which", ["lam0_min", "rho10_max", "yb0_min"])
def test_verifier3_reproduces_dual_bound(small, which):
    fr, br, T = small
    m = LKRModel2(fr, br, T)
    ver = Verifier3(fr, br, T)
    lay = m.lay
    q = np.zeros(m.nvar)
    if which == "lam0_min":
        qd = {int(lay.lam[0]): 1.0}
    elif which == "rho10_max":
        qd = {int(lay.lam[10]): -1.0, int(lay.kappa[10]): 1.0}
    else:
        qd = {int(lay.yb[0]): 1.0}
    for v, cf in qd.items():
        q[v] = cf
    val, x, z = m.solve_dual(q)
    n_eq, n_in = m._n_eq, m._n_in
    assert np.abs(m.A.T @ z + q).max() < 1e-7
    lb = ver.certify(z[:n_eq], z[n_eq:n_eq + n_in], z[n_eq + n_in:], qd, verbose=False)
    assert lb <= val + 1e-9
    assert abs(lb - val) < 1e-5, (lb, val)
