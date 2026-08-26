"""Known-answer tests: flat LCDM fits to the frozen data (FORMULATION.md §5)."""

from __future__ import annotations

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize

from . import C_KM_S
from .data import BAO, SN
from .model import whitener


def _E(z, om):
    return np.sqrt(om * (1 + z) ** 3 + (1 - om))


def dm_over_rd(z, om, h_rd):
    """D_M/r_d for flat LCDM, exact quadrature. H = 100 h E(z) km/s/Mpc, h*r_d = h_rd Mpc."""
    z = np.atleast_1d(z)
    out = np.array([quad(lambda zz: 1.0 / _E(zz, om), 0.0, zi, epsabs=0, epsrel=1e-10)[0] for zi in z])
    return C_KM_S / (100.0 * h_rd) * out


def dh_over_rd(z, om, h_rd):
    return C_KM_S / (100.0 * h_rd) / _E(np.asarray(z), om)


def chi2_bao(theta, bao: BAO, Wb):
    om, h_rd = theta
    pred = np.array([dm_over_rd(z, om, h_rd)[0] if k == "DM_over_rs" else dh_over_rd(z, om, h_rd)
                     for z, k in zip(bao.z, bao.kind)])
    r = bao.value - pred
    return float(np.sum((Wb @ r) ** 2))


def fit_bao(bao: BAO):
    Wb = whitener(bao.cov)
    res = minimize(chi2_bao, x0=[0.3, 100.0], args=(bao, Wb), method="Nelder-Mead",
                   options=dict(xatol=1e-6, fatol=1e-8, maxiter=2000))
    # crude curvature errors
    return res.x, res.fun


def fit_bao_sn(bao: BAO, sn: SN):
    """Joint fit (om, h_rd) with the SN nuisance Mp profiled analytically."""
    Wb = whitener(bao.cov)
    Wsn = whitener(sn.cov)
    ones = np.ones(len(sn.m))
    Cinv1 = Wsn.T @ (Wsn @ ones)
    denom = Cinv1 @ ones

    def chi2(theta):
        om, h_rd = theta
        D = dm_over_rd(sn.zHD, om, h_rd)
        y = sn.m - 5 * np.log10((1 + sn.zHEL) * D)
        Mp = (Cinv1 @ y) / denom
        r = y - Mp
        return chi2_bao(theta, bao, Wb) + float(np.sum((Wsn @ r) ** 2))

    res = minimize(chi2, x0=[0.3, 100.0], method="Nelder-Mead",
                   options=dict(xatol=1e-6, fatol=1e-8, maxiter=2000))
    return res.x, res.fun
