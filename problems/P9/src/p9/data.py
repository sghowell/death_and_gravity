"""Frozen-input loaders. Every file is checked against data/MANIFEST.json before use."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path(__file__).resolve().parents[2] / "data"
RAW = DATA / "raw"
CACHE = DATA / "cache"


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_manifest() -> dict:
    """Check every manifest file's sha256. Raises if any mismatch; returns the manifest."""
    man = json.loads((DATA / "MANIFEST.json").read_text())
    for entry in man["files"]:
        p = DATA / entry["path"]
        if not p.exists():
            raise FileNotFoundError(f"{p} missing; re-download from {entry['source_url']}")
        got = _sha256(p)
        if got != entry["sha256"]:
            raise ValueError(f"sha256 mismatch for {p}: {got} != {entry['sha256']}")
    return man


@dataclass(frozen=True)
class BAO:
    z: np.ndarray        # redshift per row
    value: np.ndarray    # measured D/r_d per row
    kind: tuple          # 'DM_over_rs' | 'DH_over_rs' | 'DV_over_rs'
    cov: np.ndarray      # covariance (rows in the same order)


@dataclass(frozen=True)
class SN:
    zHD: np.ndarray
    zHEL: np.ndarray
    m: np.ndarray        # m_b_corr
    cov: np.ndarray      # STAT+SYS block for the selected rows
    index: np.ndarray    # row indices into the released table


def load_desi(drop_dv: bool = True) -> BAO:
    rows = []
    for line in (RAW / "desi_gaussian_bao_ALL_GCcomb_mean.txt").read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        z, v, q = line.split()
        rows.append((float(z), float(v), q))
    cov = np.loadtxt(RAW / "desi_gaussian_bao_ALL_GCcomb_cov.txt")
    assert cov.shape == (len(rows), len(rows))
    keep = [i for i, r in enumerate(rows) if not (drop_dv and r[2] == "DV_over_rs")]
    return BAO(
        z=np.array([rows[i][0] for i in keep]),
        value=np.array([rows[i][1] for i in keep]),
        kind=tuple(rows[i][2] for i in keep),
        cov=cov[np.ix_(keep, keep)],
    )


def load_pantheon_cov() -> np.ndarray:
    """Full 1701x1701 STAT+SYS covariance (cached as .npy after first parse)."""
    CACHE.mkdir(exist_ok=True)
    npy = CACHE / "pantheon_plus_statsys.npy"
    if npy.exists():
        cov = np.load(npy)
        # the cache is derived data: re-check it against the manifest-pinned raw file's parse hash
        man = json.loads((DATA / "MANIFEST.json").read_text())
        expected = man.get("derived", {}).get("pantheon_plus_statsys_npy_sha256")
        if expected is not None and _sha256(npy) != expected:
            raise ValueError("Pantheon+ covariance cache does not match the manifest; delete data/cache and re-run")
        return cov
    with open(RAW / "Pantheon+SH0ES_STAT+SYS.cov") as f:
        n = int(f.readline())
        arr = np.loadtxt(f)
    cov = arr.reshape(n, n)
    np.save(npy, cov)
    return cov


def load_pantheon(zmin: float = 0.01, exclude_calibrators: bool = True) -> SN:
    d = pd.read_csv(RAW / "Pantheon+SH0ES.dat", sep=r"\s+")
    sel = d.zHD.values > zmin
    if exclude_calibrators:
        sel &= d.IS_CALIBRATOR.values == 0
    idx = np.where(sel)[0]
    cov = load_pantheon_cov()[np.ix_(idx, idx)]
    return SN(
        zHD=d.zHD.values[idx].astype(float),
        zHEL=d.zHEL.values[idx].astype(float),
        m=d.m_b_corr.values[idx].astype(float),
        cov=cov,
        index=idx,
    )
