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


def verify_manifest(require_sample: str | None = None) -> dict:
    """Check every manifest file's sha256. Raises if any mismatch; returns the manifest. Files marked
    "optional" (the SN-sample variants) may be absent unless they belong to require_sample."""
    man = json.loads((DATA / "MANIFEST.json").read_text())
    for entry in man["files"]:
        p = DATA / entry["path"]
        if not p.exists():
            if entry.get("optional") and entry.get("sample") != require_sample:
                continue
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
    m: np.ndarray        # m_b_corr (Pantheon+) or the released distance modulus (DES-SN5YR, Union3)
    cov: np.ndarray      # STAT+SYS block for the selected rows (for samples released as a precision matrix:
                         # its numerical inverse, used only for subsets, the LCDM fit and diagnostics)
    index: np.ndarray    # row indices into the released table
    precision: np.ndarray | None = None   # released inverse covariance (symmetrized), when that is the released
                                          # object; the certified whitening then uses it directly (model.sn_whitener)


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


# ---------------------------------------------------------------------------
# SN sample variants (FORMULATION.md §6, v2 items). Each maps to the same SN dataclass; the
# likelihood form m - 5 log10[(1+zHEL) D_M(zHD)/r_d] - M' is unchanged, M' absorbing the
# sample's magnitude/offset convention (Pantheon+: m_b_corr; DES, Union3: released mu with an
# arbitrary constant). See MANIFEST.json for the release/commit of every file.
# ---------------------------------------------------------------------------

SN_SAMPLES = ("pantheon", "dessn5yr", "union3")


def _read_snana_table(path: Path) -> pd.DataFrame:
    """SNANA key-value text table: one 'VARNAMES:' header line, then 'SN:' rows; '#' comments."""
    hdr, rows = None, []
    for line in path.read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if parts[0] == "VARNAMES:":
            hdr = parts[1:]
        elif parts[0] == "SN:":
            rows.append(parts[1:])
    df = pd.DataFrame(rows, columns=hdr)
    for c in df.columns:
        try:
            df[c] = df[c].astype(float)
        except ValueError:
            pass
    return df


def _symmetric_inverse(P: np.ndarray) -> np.ndarray:
    """Numerical cov = inv(precision), symmetrized (diagnostic / subset / LCDM-fit use only: the certified
    whitening never forms this inverse, see model.sn_whitener)."""
    P = 0.5 * (P + P.T)
    C = np.linalg.inv(P)
    return 0.5 * (C + C.T)


def load_dessn5yr_precision() -> np.ndarray:
    """DES-SN5YR (Dovekie release) STAT+SYS *inverse* covariance, 1820x1820, stored upper-triangular
    in the npz (keys 'nsn', 'cov'); unpacked exactly as in the release's DES-Dovekie-SN_Likelihood.py."""
    d = np.load(RAW / "DES-Dovekie_STAT+SYS.npz")
    n = int(d["nsn"][0])
    P = np.zeros((n, n))
    P[np.triu_indices(n)] = d["cov"]
    il = np.tril_indices(n, -1)
    P[il] = P.T[il]
    return P


def load_dessn5yr(zmin: float = 0.01) -> SN:
    """DES-SN5YR Hubble diagram (Dovekie release, 1820 SNe, zHD in [0.025, 1.14]): m := MU (bias- and
    contamination-corrected distance modulus quoted for a fixed M_0; the constant is absorbed by M'),
    cov := inv(STAT+SYS precision) restricted to zHD > zmin (no row is dropped for zmin <= 0.025)."""
    df = _read_snana_table(RAW / "DES-Dovekie_HD.csv")
    P = load_dessn5yr_precision(); P = 0.5 * (P + P.T)
    C = _symmetric_inverse(P)
    assert C.shape == (len(df), len(df)), (C.shape, len(df))
    idx = np.where(df.zHD.values > zmin)[0]
    full = len(idx) == len(df)      # the released precision is that of the full sample; a subset's is not a sub-block
    return SN(zHD=df.zHD.values[idx].astype(float), zHEL=df.zHEL.values[idx].astype(float),
              m=df.MU.values[idx].astype(float), cov=C[np.ix_(idx, idx)], index=idx, precision=P if full else None)


def _read_fits_primary_image(path: Path) -> np.ndarray:
    """Minimal reader for a FITS primary HDU holding a 2-d float64 image (no astropy dependency)."""
    raw = path.read_bytes()
    cards, pos = {}, 0
    while True:
        card = raw[pos:pos + 80].decode("ascii"); pos += 80
        key = card[:8].strip()
        if key == "END":
            break
        if card[8:10] == "= ":
            cards[key] = card[10:].split("/")[0].strip()
    if cards.get("BITPIX") != "-64" or cards.get("NAXIS") != "2":
        raise ValueError(f"unexpected FITS layout {cards}")
    n1, n2 = int(cards["NAXIS1"]), int(cards["NAXIS2"])
    start = ((pos + 2879) // 2880) * 2880
    return np.frombuffer(raw[start:start + 8 * n1 * n2], dtype=">f8").reshape(n2, n1).astype(float)


def load_union3() -> SN:
    """Union3 + UNITY1.5 compressed Hubble diagram (Rubin et al. 2023): 22 binned distance moduli.
    File layout (release README): first row = redshifts, first column = mu, rest = inverse covariance.
    The mu's carry an arbitrary additive constant (degenerate with the UNITY script-M; the release
    subtracts a per-sample median), absorbed by M'. Bins are in the CMB frame; zHEL := zHD := z_bin, so
    D_L = (1+z) D_M(z) at the bin redshift (the (1+zHEL) factor is the bin's own redshift)."""
    a = _read_fits_primary_image(RAW / "mu_mat_union3_cosmo=2_mu.fits")
    z = a[0, 1:].copy(); mu = a[1:, 0].copy(); icov = a[1:, 1:]; icov = 0.5 * (icov + icov.T)
    return SN(zHD=z, zHEL=z.copy(), m=mu, cov=_symmetric_inverse(icov), index=np.arange(len(z)), precision=icov)


def load_sn(name: str = "pantheon", zmin: float = 0.01) -> SN:
    """Dispatcher over SN_SAMPLES (zmin applies to the per-SN samples; Union3 bins start at z = 0.05)."""
    if name == "pantheon":
        return load_pantheon(zmin=zmin)
    if name == "dessn5yr":
        return load_dessn5yr(zmin=zmin)
    if name == "union3":
        return load_union3()
    raise ValueError(f"unknown SN sample {name!r}; choose from {SN_SAMPLES}")
