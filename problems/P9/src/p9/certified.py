"""Certified version of the bound computation (FORMULATION.md §4).

Every optimization whose value enters the final claim is certified by a conic dual vector
checked in ball arithmetic (verify.Verifier): the node-bound tightening passes (on the SN subset
model, a valid outer set) and the final bound (full sample). Node bounds passed between passes are
the rigorous (outward-rounded) values. T is the rigorous upper enclosure of chi2 at the reference
point plus Delta. Certificates are stored under results/certificates/<tag>/ as .npz files.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from flint import arb

from . import C_KM_S
from .certify import solve_with_dual
from .classmin import in_class, minimize_chi2_over_class
from .model import Frozen
from .socp2 import KappaModel, NodeBounds
from .verify import Verifier, _endpoint, rigorous_chi2

RESULTS = Path(__file__).resolve().parents[2] / "results"


def _save_cert(path: Path, cert, extra: dict):
    np.savez_compressed(path, q_u=cert.q_u, z_eq=cert.z_eq, z_in=cert.z_in, z_soc=cert.z_soc,
                        primal_obj=cert.primal_obj, dual_obj=cert.dual_obj, **extra)


def certified_node_bounds(fr: Frozen, nb: NodeBounds, T: float, tangent_u: np.ndarray,
                          cert_dir: Path, pass_id: int, verbose: bool = True) -> NodeBounds:
    """Rigorous min/max of u_i and D_i over the relaxed set built from nb (subset model)."""
    m = KappaModel(fr, nb, T, tangent_u=tangent_u)
    ver = Verifier(fr, nb, T, tangent_u)
    N = fr.spec.n_seg
    u_lo = np.empty(N + 1); u_hi = np.empty(N + 1); dm_lo = np.zeros(N + 1); dm_hi = np.zeros(N + 1)
    t0 = time.time()
    for i in range(N + 1):
        e = np.zeros(N + 1); e[i] = 1.0
        for name, q_u in (("u_lo", e), ("u_hi", -e), ("dm_lo", fr.A_nodes[i]), ("dm_hi", -fr.A_nodes[i])):
            if i == 0 and name.startswith("dm"):
                continue
            cert = solve_with_dual(m, q_u)
            lb = ver.certify(cert.z_eq, cert.z_in, cert.z_soc, q_u, verbose=False)
            _save_cert(cert_dir / f"pass{pass_id}_node{i}_{name}.npz", cert, dict(rigorous=lb, pass_id=pass_id, node=i))
            if name == "u_lo":
                u_lo[i] = lb
            elif name == "u_hi":
                u_hi[i] = -lb
            elif name == "dm_lo":
                dm_lo[i] = lb
            else:
                dm_hi[i] = -lb
        if verbose and i % 15 == 0:
            print(f"    [{time.time()-t0:5.0f}s] node {i:3d}: u in [{u_lo[i]:.5f},{u_hi[i]:.5f}] "
                  f"D in [{dm_lo[i]:.5f},{dm_hi[i]:.5f}]", flush=True)
    # intersect with the previous (also rigorous) bounds
    return NodeBounds(np.maximum(u_lo, nb.u_lo), np.minimum(u_hi, nb.u_hi),
                      np.maximum(dm_lo, nb.dm_lo), np.minimum(dm_hi, nb.dm_hi))


def certified_bound(fr_full: Frozen, fr_sub: Frozen, Delta: float, u_start: np.ndarray, tag: str,
                    n_passes: int = 6, tol_rel: float = 1e-4, verbose: bool = True) -> dict:
    sp = fr_full.spec
    cert_dir = RESULTS / "certificates" / tag
    cert_dir.mkdir(parents=True, exist_ok=True)
    # reference point: polished class minimizer; T from a rigorous enclosure of its chi2
    u_ref, Mp_ref, chi2_ref = minimize_chi2_over_class(fr_full, u_start)
    assert in_class(fr_full, u_ref)
    chi2_ball = rigorous_chi2(fr_full, u_ref, Mp_ref)
    chi2_up = _endpoint(chi2_ball, +1)
    T = chi2_up + Delta
    if verbose:
        print(f"  reference: chi2 (float) {chi2_ref:.6f}, rigorous {chi2_ball.str(12)} -> T = {T:.6f}", flush=True)
    nb = NodeBounds.a_priori(fr_full)          # definition of the class: no certificate needed
    N = sp.n_seg
    e0 = np.zeros(N + 1); e0[0] = 1.0
    history = []
    last = None
    lb = None
    for it in range(n_passes):
        # (a) improve the reference point (any class point is valid; only T changes)
        m_full = KappaModel(fr_full, nb, T, tangent_u=u_ref)
        _, u_c, _ = m_full.min_chi2()
        u_c, Mp_c, chi2_c = minimize_chi2_over_class(fr_full, u_c)
        if chi2_c < chi2_ref - 1e-9:
            u_ref, Mp_ref, chi2_ref = u_c, Mp_c, chi2_c
            chi2_ball = rigorous_chi2(fr_full, u_ref, Mp_ref)
            T = _endpoint(chi2_ball, +1) + Delta
            m_full = KappaModel(fr_full, nb, T, tangent_u=u_ref)
        # (b) certified bound with the full sample
        ver_full = Verifier(fr_full, nb, T, u_ref)
        cert = solve_with_dual(m_full, e0)
        lb = ver_full.certify(cert.z_eq, cert.z_in, cert.z_soc, e0, verbose=False)
        _save_cert(cert_dir / f"pass{it}_final_u0.npz", cert, dict(rigorous=lb, pass_id=it, T=T,
                                                                  u_ref=u_ref, Mp_ref=Mp_ref))
        H0 = C_KM_S / (sp.r_lo * lb)
        history.append(dict(pass_=it, T=T, chi2_ref=chi2_ref, u0_min=lb, solver_u0=cert.primal_obj,
                            H0_max=H0, width=nb.width()))
        if verbose:
            print(f"  pass {it}: chi2_ref={chi2_ref:.4f} T={T:.4f} rigorous u0_min={lb:.6f} "
                  f"(solver {cert.primal_obj:.6f}) H0_max={H0:.4f} width={nb.width():.4f}", flush=True)
        if last is not None and abs(lb - last) <= tol_rel * last:
            break
        last = lb
        # (c) certified node bounds on the subset model
        nb = certified_node_bounds(fr_sub, nb, T, u_ref, cert_dir, it, verbose=verbose)
    result = dict(tag=tag, L=sp.L, Delta=Delta, T=T, chi2_ref=chi2_ref, u0_min=lb,
                  H0_max=C_KM_S / (sp.r_lo * lb), r_lo=sp.r_lo, H0_ref=C_KM_S / (147.09 * u_ref[0]),
                  u_ref=u_ref.tolist(), Mp_ref=Mp_ref, history=history,
                  node_bounds=dict(u_lo=nb.u_lo.tolist(), u_hi=nb.u_hi.tolist(),
                                   dm_lo=nb.dm_lo.tolist(), dm_hi=nb.dm_hi.tolist()),
                  grid_x=sp.x.tolist(), sn_subset_index=fr_sub.sn.index.tolist())
    (cert_dir / "result.json").write_text(json.dumps(result, indent=1))
    return result
