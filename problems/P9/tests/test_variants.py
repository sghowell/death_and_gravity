"""FORMULATION §6.1 options: SN-sample loaders (whitening from the released precision), the BGS D_V row in
the shared row builder (verifier reproduces solver duals; every row holds at random class members), and
the r_d box / variant plumbing. Small subsets only."""
import numpy as np
import pytest
from flint import arb
from p9 import C_KM_S
from p9.bound import sn_subset
from p9.classmin import in_class
from p9.data import RAW, load_desi, load_sn
from p9.geometry import lcdm_u_nodes
from p9.lcdm import fit_bao
from p9.lkr import Brackets3, initial_brackets3
from p9.lkr2 import LKRModel2
from p9.lkr_rows import FloatArith, build, layout_for, obbt_objectives
from p9.model import ClassSpec, Frozen, sn_whitener
from p9.variants import RD_BOX_PLANCK, Variant, parse_rd_box, variant_tag
from p9.verify import rigorous_chi2
from p9.verify3 import Verifier3

FILES = {"dessn5yr": ["DES-Dovekie_HD.csv", "DES-Dovekie_STAT+SYS.npz"], "union3": ["mu_mat_union3_cosmo=2_mu.fits"]}


def _need(name):
    if not all((RAW / f).exists() for f in FILES[name]):
        pytest.skip(f"{name} raw files not downloaded (see MANIFEST.json)")


# ---------------------------------------------------------------- SN samples
@pytest.mark.parametrize("name,n_expected", [("dessn5yr", 1820), ("union3", 22)])
def test_precision_released_samples_whiten_without_inverse(name, n_expected):
    _need(name)
    sn = load_sn(name); n = len(sn.m)
    assert n == n_expected and sn.precision is not None and sn.cov.shape == (n, n)
    W = sn_whitener(sn)
    assert np.allclose(W, np.triu(W))                                  # W = L^T, P = L L^T
    assert np.abs(W.T @ W - sn.precision).max() <= 1e-12 * np.abs(sn.precision).max()
    assert np.abs(W @ sn.cov @ W.T - np.eye(n)).max() < 1e-10          # C = inv(P) only enters this diagnostic
    assert np.all(sn.zHD > 0.01) and sn.zHD.max() < 2.5 and np.all(sn.zHEL > 0)
    # a Hubble diagram: mu rises with z (weighted by 1/sigma^2; some DES SNe have BEAMS-inflated sigmas of ~20 mag)
    w = 1.0 / np.diag(sn.cov); lo, hi = sn.zHD < np.median(sn.zHD), sn.zHD >= np.median(sn.zHD)
    assert np.average(sn.m[hi], weights=w[hi]) > np.average(sn.m[lo], weights=w[lo]) + 0.5


def test_pantheon_path_unchanged():
    sn = load_sn("pantheon")
    assert sn.precision is None and len(sn.m) == 1580


def test_union3_offset_is_absorbed_by_Mp():
    """The Union3 mu's carry an arbitrary constant; the profiled M' of the LCDM point must stay inside MP_BOX."""
    _need("union3")
    from p9.socp2 import MP_BOX
    bao = load_desi(); sn = load_sn("union3"); spec = ClassSpec(L=1.5, grid_kind="geometric"); fr = Frozen(bao, sn, spec)
    Mp = fr.best_Mp(lcdm_u_nodes(spec.x, 0.30, 101.0))
    assert MP_BOX[0] + 1 < Mp < MP_BOX[1] - 1


# ---------------------------------------------------------------- D_V row
def test_bao_only_lcdm_with_dv_reproduces_desi_dr2():
    """Known-answer test 1 (FORMULATION §5) with the full 13-row vector, D_V = (z D_M^2 D_H)^{1/3}."""
    (om, hrd), _ = fit_bao(load_desi(drop_dv=False))
    assert abs(om - 0.2975) < 0.002 and abs(hrd - 101.54) < 0.05


@pytest.fixture(scope="module")
def small_dv():
    bao = load_desi(drop_dv=False); sub = sn_subset(load_sn("pantheon"), 40)
    spec = ClassSpec(L=1.5, grid_kind="geometric", use_dv=True); fr = Frozen(bao, sub, spec)
    u = lcdm_u_nodes(spec.x, 0.30, 101.0)
    T = fr.chi2(u, fr.best_Mp(u)) + 4.0
    return fr, initial_brackets3(fr), T, u


def test_layout_is_baseline_when_dv_off(small_dv):
    fr, br, T, u = small_dv
    fr0 = Frozen(load_desi(), fr.sn, ClassSpec(L=1.5, grid_kind="geometric"))
    lay0, lay1 = layout_for(fr0), layout_for(fr)
    assert lay0.idx_dv == [] and len(lay0.yb) == 6 and len(lay0.yH) == 0 and len(lay0.yV) == 0
    assert lay1.idx_dv == [0] and len(lay1.yb) == 7 and len(lay1.yH) == 1 and lay1.dvnodes == [42, 43]
    assert lay1.nvar == lay0.nvar + 5          # y_b, y_H, y_V, P_V, w_b for the extra row
    assert len(obbt_objectives(lay1)) == len(obbt_objectives(lay0)) + 2 * (1 + 1 + 2)   # yb, yv, two lambda nodes


def test_rigorous_chi2_with_dv_matches_float(small_dv):
    fr, br, T, u = small_dv
    Mp = fr.best_Mp(u); ball = rigorous_chi2(fr, u, Mp)
    assert abs(float(ball.mid().str(20, radius=False)) - fr.chi2(u, Mp)) < 1e-8
    assert fr.dv_rows == [0] and abs(fr.bao_pred(u)[0] - np.cbrt(fr.bao.z[0] * (fr.AM_bao[0] @ u) ** 2 * (fr.AH_bao[0] @ u))) < 1e-12


@pytest.mark.parametrize("which", ["lam0_min", "yv0_min", "yv0_max", "rho10_max", "yb6_min"])
def test_verifier3_reproduces_dual_bound_with_dv(small_dv, which):
    fr, br, T, u = small_dv
    m = LKRModel2(fr, br, T); ver = Verifier3(fr, br, T); lay = m.lay
    qd = {"lam0_min": {int(lay.lam[0]): 1.0}, "yv0_min": {int(lay.yV[0]): 1.0}, "yv0_max": {int(lay.yV[0]): -1.0},
          "rho10_max": {int(lay.lam[10]): -1.0, int(lay.kappa[10]): 1.0}, "yb6_min": {int(lay.yb[6]): 1.0}}[which]
    q = np.zeros(m.nvar)
    for v, cf in qd.items():
        q[v] = cf
    val, x, z = m.solve_dual(q)
    n_eq, n_in = m._n_eq, m._n_in
    assert np.abs(m.A.T @ z + q).max() < 1e-7
    lb = ver.certify(z[:n_eq], z[n_eq:n_eq + n_in], z[n_eq + n_in:], qd, verbose=False)
    assert lb <= val + 1e-9 and abs(lb - val) < 1e-5, (lb, val)


def _assignment(fr, lay, u, Mp):
    """Every relaxation variable evaluated exactly (floats) at a class member u."""
    x = fr.spec.x; N = lay.N
    lam = np.log10(u); D = fr.A_nodes @ u
    kappa = np.concatenate([[lam[0]], np.log10(D[1:] / np.expm1(x[1:]))])
    v = np.zeros(lay.nvar)
    v[lay.lam] = lam; v[lay.kappa] = kappa; v[lay.Mp] = Mp
    v[lay.ell] = np.log10(fr.A_sn @ u)
    v[lay.Ed[0]] = D[1] / u[0]; v[lay.Es[0]] = u[1] / u[0]
    for i in range(1, N):
        v[lay.Ed[i]] = D[i + 1] / D[i]; v[lay.Er[i]] = u[i] / D[i]; v[lay.Es[i]] = u[i + 1] / D[i]
    pred = fr.bao_pred(u)
    v[lay.P] = pred
    for p, r in enumerate(lay.idx_ym):
        v[lay.yb[p]] = np.log10(fr.AM_bao[r] @ u)
    for q, r in enumerate(lay.idx_dv):
        v[lay.yH[q]] = np.log10(fr.AH_bao[r] @ u); v[lay.yV[q]] = np.log10(pred[r])
    for node in lay.enodes:
        v[lay.U[lay.epos[node]]] = u[node]
    v[lay.wb] = fr.Wb @ (fr.bao.value - pred)
    v[lay.ws] = fr.Wsn @ (fr.sn.m - fr.sn_offset - 5.0 * v[lay.ell] - Mp)
    return v


def test_relaxation_rows_hold_at_random_class_members(small_dv):
    """Every eq/le row of the shared builder (D_V rows included) and every absorption bound must hold at
    random members of C(G, L) with the class-box brackets: D_V computed exactly from u."""
    fr, br, T, _ = small_dv
    spec = fr.spec; N = spec.n_seg; lh = np.log10(1 + spec.L * spec.hs)
    B = build(fr, br, None, FloatArith()); lay = B.layout
    rng = np.random.default_rng(1)
    worst_eq = worst_le = 0.0
    for trial in range(60):
        lam0 = np.log10(C_KM_S / (rng.uniform(40, 120) * 147.09))
        steps = rng.uniform(-1, 1, N) * lh * rng.choice([1.0, 0.999999, 0.3])
        u = 10 ** (lam0 + np.concatenate([[0.0], np.cumsum(steps)]))
        u = np.clip(u, spec.u_box[0] * 1.001, spec.u_box[1] * 0.999)
        if not in_class(fr, u):
            continue
        v = _assignment(fr, lay, u, fr.best_Mp(u))
        for coeffs, rhs in B.eq_rows:
            worst_eq = max(worst_eq, abs(sum(cf * v[i] for i, cf in coeffs) - rhs) / max(1.0, abs(rhs)))
        for coeffs, rhs in B.le_rows:
            worst_le = max(worst_le, (sum(cf * v[i] for i, cf in coeffs) - rhs) / max(1.0, abs(rhs)))
        for k in range(lay.nvar):
            assert B.var_lo[k] - 1e-9 <= v[k] <= B.var_hi[k] + 1e-9, (k, v[k], B.var_lo[k], B.var_hi[k])
    assert worst_eq < 1e-9 and worst_le < 1e-9, (worst_eq, worst_le)


def test_tighten_with_dv_keeps_lcdm_point(small_dv):
    fr, br, T, u = small_dv
    m = LKRModel2(fr, br, T)
    br2 = m.tighten(rho_nodes=range(1, 6), lam_nodes=[0] + m.lay.dvnodes)
    yV = np.log10(fr.bao_pred(u)[fr.dv_rows[0]])
    assert len(br2.yv_lo) == 1 and br2.yv_lo[0] - 1e-9 <= yV <= br2.yv_hi[0] + 1e-9
    assert br2.yv_hi[0] - br2.yv_lo[0] < br.yv_hi[0] - br.yv_lo[0]
    rt = Brackets3.from_dict(br2.to_dict())
    assert all(np.array_equal(getattr(rt, k), getattr(br2, k)) for k in Brackets3.KEYS)
    # apply_bound intersects (never loosens) and rounds an Arb shift outward
    b = br2.copy(); old = b.rho_lo[3]
    b.apply_bound("rho", 3, "lo", old + 1.25, arb(0.25)); assert old + 1.0 - 1e-12 < b.rho_lo[3] <= old + 1.0
    b.apply_bound("rho", 3, "lo", old, arb(0.0)); assert b.rho_lo[3] > old + 1.0 - 1e-12          # weaker bound: no-op
    old = b.yv_hi[0]; b.apply_bound("yv", 0, "hi", -(old - 1e-3), 0.0); assert abs(b.yv_hi[0] - (old - 1e-3)) < 1e-12


# ---------------------------------------------------------------- r_d box and variants
def test_variant_tags_and_rd_box():
    assert variant_tag("pantheon", False, *RD_BOX_PLANCK) == "" and Variant.from_state({"r_lo": 146.57}).tag == ""
    assert Variant("union3", True, *RD_BOX_PLANCK).tag == "_union3_dv"
    lo, hi = parse_rd_box(["bbn"])
    assert 145 < lo < 147.09 < hi < 151 and parse_rd_box(["planck"]) == RD_BOX_PLANCK
    assert parse_rd_box(["145.5", "150.5"]) == (145.5, 150.5)
    assert Variant("pantheon", False, lo, hi).tag == f"_rd{lo:g}-{hi:g}"
    # the class box follows the r_d box; the certified bound rescales exactly as 1/r_d
    s0 = ClassSpec(L=1.5, grid_kind="geometric"); s1 = ClassSpec(L=1.5, grid_kind="geometric", r_lo=lo, r_hi=hi)
    assert s1.u_box[0] < s0.u_box[0] and s1.u_box[1] > s0.u_box[1]
    u0_min = 0.5
    assert abs(C_KM_S / (s0.r_lo * u0_min) * s0.r_lo / 147.09 - C_KM_S / (147.09 * u0_min)) < 1e-9
