"""Variant selection shared by the drivers (FORMULATION.md §6.1): SN sample, BGS D_V row, r_d box.

The baseline (Pantheon+, no D_V row, Planck ±2σ r_d box) has the empty variant tag, so the certificate
directories and result files of the baseline chains are unchanged; every variant gets a suffix
'_<sn>' + '_dv' + '_rd<lo>-<hi>' (only the non-default parts).
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass

from .data import DATA, SN_SAMPLES, load_desi, load_sn, verify_manifest
from .model import ClassSpec, Frozen

RD_BOX_PLANCK = (ClassSpec.r_lo, ClassSpec.r_hi)      # (146.57, 147.61): Planck 2018 r_drag ± 2σ


def rd_box_named(name: str) -> tuple[float, float]:
    """Named r_d boxes: 'planck' (the ClassSpec default) or 'bbn' (MANIFEST scalar_inputs / r_drag_BBN)."""
    if name == "planck":
        return RD_BOX_PLANCK
    if name == "bbn":
        man = json.loads((DATA / "MANIFEST.json").read_text())
        e = next(s for s in man["scalar_inputs"] if s["name"] == "r_drag_BBN")
        return float(e["box_Mpc"][0]), float(e["box_Mpc"][1])
    raise ValueError(f"unknown r_d box {name!r}")


def parse_rd_box(v) -> tuple[float, float]:
    if v is None:
        return RD_BOX_PLANCK
    if len(v) == 1:
        return rd_box_named(v[0])
    if len(v) == 2:
        lo, hi = float(v[0]), float(v[1])
        if not (0 < lo <= hi):
            raise ValueError(f"bad r_d box {v}")
        return lo, hi
    raise ValueError(f"--rd_box takes 'planck', 'bbn' or two numbers, got {v}")


def add_variant_args(ap: argparse.ArgumentParser) -> None:
    ap.add_argument("--sn", choices=SN_SAMPLES, default="pantheon", help="SN sample (default pantheon)")
    ap.add_argument("--dv", action="store_true", help="include the DESI BGS D_V/r_d row (ClassSpec.use_dv)")
    ap.add_argument("--rd_box", nargs="+", default=None, metavar="LO_HI",
                    help="r_d box [Mpc]: two numbers, or 'planck' (default: Planck 2018 ±2σ) / 'bbn' (MANIFEST r_drag_BBN)")


def variant_tag(sn: str, use_dv: bool, r_lo: float, r_hi: float) -> str:
    tag = "" if sn == "pantheon" else f"_{sn}"
    if use_dv:
        tag += "_dv"
    if (float(r_lo), float(r_hi)) != RD_BOX_PLANCK:
        tag += f"_rd{r_lo:g}-{r_hi:g}"
    return tag


@dataclass(frozen=True)
class Variant:
    sn: str = "pantheon"
    use_dv: bool = False
    r_lo: float = RD_BOX_PLANCK[0]
    r_hi: float = RD_BOX_PLANCK[1]

    @staticmethod
    def from_args(a) -> Variant:
        lo, hi = parse_rd_box(a.rd_box)
        return Variant(a.sn, bool(a.dv), lo, hi)

    @staticmethod
    def from_state(s: dict) -> Variant:
        """From a chain's state.json (baseline chains predate the variant keys: defaults apply)."""
        return Variant(s.get("sn", "pantheon"), bool(s.get("use_dv", False)),
                       float(s.get("r_lo", RD_BOX_PLANCK[0])), float(s.get("r_hi", RD_BOX_PLANCK[1])))

    @property
    def tag(self) -> str:
        return variant_tag(self.sn, self.use_dv, self.r_lo, self.r_hi)

    def as_dict(self) -> dict:
        return dict(sn=self.sn, use_dv=self.use_dv, r_lo=self.r_lo, r_hi=self.r_hi)

    def describe(self) -> str:
        return f"sn={self.sn} dv={self.use_dv} r_d in [{self.r_lo:g}, {self.r_hi:g}] Mpc"

    def frozen(self, L: float, refine: int, grid_kind: str = "geometric"):
        """(bao, sn, spec, Frozen) for this variant; the manifest is verified (the variant's files required)."""
        verify_manifest(require_sample=self.sn)
        bao = load_desi(drop_dv=not self.use_dv)
        sn = load_sn(self.sn)
        spec = ClassSpec(L=L, grid_kind=grid_kind, refine=refine, use_dv=self.use_dv, r_lo=self.r_lo, r_hi=self.r_hi)
        return bao, sn, spec, Frozen(bao, sn, spec)
