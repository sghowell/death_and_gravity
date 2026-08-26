"""Ball-arithmetic verifier for LKR dual certificates, consuming the SAME row builder as the solver
(lkr_rows.build with ArbArith). Weak duality with residual absorption, as in verify.py:

    q^T y >= -b^T z - sum_k |rho_k| Y_k,   rho = A^T z + q,  |y_k| <= Y_k on the feasible set,
provided z_le >= 0 and z_soc in the second-order cone. Row order: eq = [BAO rows | SN dense block |
remaining eq rows], le = build order, soc = (sqrt T, w_b, w_s).
"""

from __future__ import annotations

import numpy as np
from flint import arb, arb_mat

from .lkr import Brackets3
from .lkr_rows import ArbArith, build
from .model import Frozen
from .verify import _endpoint


class Verifier3:
    def __init__(self, fr: Frozen, br: Brackets3, T: float):
        self.fr = fr
        self.ar = ArbArith()
        self.B = build(fr, br, T, self.ar)
        self.lay = self.B.layout
        self.Wsn = arb_mat([[arb(float(v)) for v in row] for row in fr.Wsn])
        self.sn_rhs = arb_mat([[v] for v in self.B.sn_rhs])

    def certify(self, z_eq: np.ndarray, z_le: np.ndarray, z_soc: np.ndarray, q: dict, verbose: bool = True) -> float:
        """q: {var_index: float coefficient} of the objective. Returns a rigorous lower bound on q.y over F_rel."""
        B, lay = self.B, self.lay
        nvar = lay.nvar; n, nbao = lay.n, lay.nbao
        ar = self.ar
        if np.any(z_le < 0):
            raise ValueError("dual inequality multipliers must be nonnegative")
        z_soc = z_soc.copy()
        nrm = arb(0)
        for v in z_soc[1:]:
            nrm += arb(float(v)) ** 2
        nrm = nrm.sqrt()
        if not (arb(float(z_soc[0])) >= nrm):
            z_soc[0] = _endpoint(nrm, +1)
        rho = [arb(0) for _ in range(nvar)]
        val = arb(0)
        # ---- equalities ----
        n_bao_rows = nbao
        eq_rows = B.eq_rows
        zpos = 0
        for r in range(n_bao_rows):
            coeffs, rhs = eq_rows[r]; zr = arb(float(z_eq[zpos])); zpos += 1
            for v, cf in coeffs:
                rho[v] += cf * zr
            val -= rhs * zr
        # SN dense block: rows w_s + 5 Wsn ell + (Wsn 1) Mp = Wsn (m - off)
        zs = arb_mat([[arb(float(v))] for v in z_eq[zpos:zpos + n]]); zpos += n
        WsnT_zs = self.Wsn.transpose() * zs
        ones = arb_mat([[arb(1)]] * n)
        Wsn1 = self.Wsn * ones
        rhs_sn = self.Wsn * self.sn_rhs
        for j in range(n):
            rho[lay.ell[j]] += 5 * WsnT_zs[j, 0]
            rho[lay.ws[j]] += zs[j, 0]
            rho[lay.Mp] += Wsn1[j, 0] * zs[j, 0]
            val -= rhs_sn[j, 0] * zs[j, 0]
        for r in range(n_bao_rows, len(eq_rows)):
            coeffs, rhs = eq_rows[r]; zr = arb(float(z_eq[zpos])); zpos += 1
            for v, cf in coeffs:
                rho[v] += cf * zr
            val -= rhs * zr
        assert zpos == len(z_eq), (zpos, len(z_eq))
        # ---- inequalities ----
        assert len(z_le) == len(B.le_rows), (len(z_le), len(B.le_rows))
        for (coeffs, rhs), zv in zip(B.le_rows, z_le):
            if zv == 0.0:
                continue
            zr = arb(float(zv))
            for v, cf in coeffs:
                rho[v] += cf * zr
            val -= rhs * zr
        # ---- SOC ----
        z0 = arb(float(z_soc[0]))
        val -= B.soc_T * z0
        for r in range(nbao):
            rho[lay.wb[r]] += -arb(float(z_soc[1 + r]))
        for j in range(n):
            rho[lay.ws[j]] += -arb(float(z_soc[1 + nbao + j]))
        # ---- objective ----
        for v, cf in q.items():
            rho[v] += arb(float(cf))
        # ---- absorption ----
        loss = arb(0)
        for k in range(nvar):
            Y = ar.absmax(B.var_lo[k], B.var_hi[k])
            loss += abs(rho[k]) * Y
        bound = val - loss
        lb = _endpoint(bound, -1)
        self._val, self._loss = val, loss
        if verbose:
            print(f"    certificate: -b^T z = {val.str(12)}  loss = {loss.str(5)}  => rigorous lower bound {lb:.8f}", flush=True)
        return lb
