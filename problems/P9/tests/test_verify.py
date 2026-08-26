"""The rigorous verifier must reproduce the solver's dual bound (up to the tiny residual loss)
and must be consistent with the model's constraint matrix."""
import numpy as np
import pytest

from p9.bound import sn_subset
from p9.certify import solve_with_dual
from p9.data import load_desi, load_pantheon
from p9.geometry import lcdm_u_nodes
from p9.model import ClassSpec, Frozen
from p9.socp2 import KappaModel, NodeBounds
from p9.verify import Verifier, rigorous_chi2


@pytest.fixture(scope="module")
def small():
    bao = load_desi(); sn = load_pantheon(); sub = sn_subset(sn, 40)
    spec = ClassSpec(L=1.5, grid_kind="geometric")
    fr = Frozen(bao, sub, spec)
    u = lcdm_u_nodes(spec.x, 0.30, 101.0)
    Mp = fr.best_Mp(u)
    T = fr.chi2(u, Mp) + 4.0
    return fr, u, T


def test_rigorous_chi2_matches_float(small):
    fr, u, T = small
    Mp = fr.best_Mp(u)
    ball = rigorous_chi2(fr, u, Mp)
    mid = float(ball.mid().str(20, radius=False))
    assert abs(mid - fr.chi2(u, Mp)) < 1e-8


@pytest.mark.parametrize("which", ["u0_min", "u0_max", "D40_min"])
def test_certificate_reproduces_solver_bound(small, which):
    fr, u, T = small
    nb = NodeBounds.a_priori(fr)
    m = KappaModel(fr, nb, T, tangent_u=u)
    ver = Verifier(fr, nb, T, u)
    N = fr.spec.n_seg
    e0 = np.zeros(N + 1); e0[0] = 1.0
    q = {"u0_min": e0, "u0_max": -e0, "D40_min": fr.A_nodes[40]}[which]
    cert = solve_with_dual(m, q)
    z = np.concatenate([cert.z_eq, cert.z_in, cert.z_soc])
    qq = np.zeros(m.nvar); qq[m.idx["u"]] = q
    assert np.abs(m.A.T @ z + qq).max() < 1e-8          # solver dual is (nearly) feasible
    lb = ver.certify(cert.z_eq, cert.z_in, cert.z_soc, q, verbose=False)
    assert lb <= cert.primal_obj + 1e-9                   # rigorous bound never exceeds the primal
    assert abs(lb - cert.primal_obj) < 1e-6              # and is tight
