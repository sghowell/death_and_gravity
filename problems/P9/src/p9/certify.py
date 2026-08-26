"""Certificate generation: solve a KappaModel instance, extract the conic dual, verify rigorously."""

from __future__ import annotations

from dataclasses import dataclass

import clarabel
import numpy as np

from .socp2 import KappaModel
from .verify import Verifier


@dataclass
class DualCertificate:
    q_u: np.ndarray
    z_eq: np.ndarray
    z_in: np.ndarray
    z_soc: np.ndarray
    primal_obj: float
    dual_obj: float


def solve_with_dual(m: KappaModel, q_u: np.ndarray, q_k: np.ndarray | None = None) -> DualCertificate:
    q = np.zeros(m.nvar); q[m.idx["u"]] = q_u
    if q_k is not None:
        q[m.idx["kappa"]] = q_k
    solver = clarabel.DefaultSolver(m.P0, q, m.A, m.b, m.cones, m.settings)
    sol = solver.solve()
    if "Solved" not in str(sol.status):
        raise RuntimeError(str(sol.status))
    z = np.asarray(sol.z)
    n_eq, n_in = m._n_eq, m._n_in
    z_eq, z_in, z_soc = z[:n_eq], z[n_eq:n_eq + n_in], z[n_eq + n_in:]
    return DualCertificate(q_u=q_u, z_eq=z_eq, z_in=z_in, z_soc=z_soc,
                           primal_obj=float(sol.obj_val), dual_obj=float(-m.b @ z))


def certify(m: KappaModel, ver: Verifier, q_u: np.ndarray, verbose: bool = True) -> tuple[float, DualCertificate]:
    """Rigorous lower bound on q_u . u over the relaxed set, with the certificate that proves it."""
    cert = solve_with_dual(m, q_u)
    lb = ver.certify(cert.z_eq, cert.z_in, np.maximum(cert.z_soc, 0) if len(cert.z_soc) else cert.z_soc,
                     q_u, verbose=verbose)
    if verbose:
        print(f"    solver primal {cert.primal_obj:.8f}  dual {cert.dual_obj:.8f}  rigorous {lb:.8f}")
    return lb, cert
