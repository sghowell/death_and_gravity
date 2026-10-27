"""S6 part (iii): pole-aware band counts for l <= 19 (pre-registration O5), KAT boxes, controls, replays.

Band box B_l = [-a'_l, a'_l] x [-b', -b0] (exact dyadics): a'_l from part (ii') (certificates/beta/part2),
b' = e kappa_- rounded up, b0 = kappa_-/2 rounded down (beta_common.band_bounds).  S3's disc machinery
(part3.plan_edge, part3.segment_with_split; pole clearance on every disc box) gives the argument
increment of W_l along the right half H of the boundary, (0, -b') -> (a', -b') -> (a', -b0) -> (0, -b0);
by the symmetry W(-conj w) = conj W(w) the winding is (sum over H)/pi.  W_l is meromorphic in B_l with
at most one pole, at omega0 = -i kappa_+, of order <= 1, so winding = Z - P with P in {0, 1}.
Expected: winding 0 at l = 0 (Z = 1: omega_*; P = 1), winding -1 for 1 <= l <= 19 (P = 1, Z = 0).

Commands (PYTHONPATH=src):
  python -m p2.beta_part3 run [workers] [retry]      (detached via launch; resumable JSONL in results/beta/part3)
  python -m p2.beta_part3 launch [workers] [retry]
  python -m p2.beta_part3 assemble                    (certificates/beta/part3/band_PB_l{l}.json)
  python -m p2.beta_part3 kat [workers]               (KAT boxes B2-B4, B6; NB1, NB2; left halves; B1)
  python -m p2.beta_part3 replay [workers]            (5% replay of the band segments)
"""
from __future__ import annotations

import itertools
import json
import math
import os
import platform
import socket
import subprocess
import sys
import time
from pathlib import Path

from flint import acb, arb, fmpq

from .arbutil import Timer, precision
from .beta_common import CERT, EDGE, L_SMALL, RES, band_bounds, commit_tag, write_json
from .params import Point

RES3 = RES / "part3"
CERT3 = CERT / "part3"
SRC = Path(__file__).resolve().parents[1]
A_FIXED = CERT / "part2" / "A_fixed_PB.json"


def a_prime(l: int) -> fmpq:
    rows = {r["l"]: r for r in json.loads(A_FIXED.read_text())["rows"]}
    q = rows[l]["a_prime"]
    n, d = (int(x) for x in q.split("/")) if "/" in q else (int(q), 1)
    return fmpq(n, d)


def corners(pt: Point, l: int, bottom: fmpq | None = None):
    bp, b0 = band_bounds(pt)
    if bottom is not None:
        bp = bottom
    a = a_prime(l)
    return [(fmpq(0), -bp), (a, -bp), (a, -b0), (fmpq(0), -b0)], a, bp, b0


def plan_path(pt: Point, cs: list, edge0: int = 0) -> list:
    from . import part3
    km = float(pt.km.mid())
    segs = []
    for k, (z0, z1) in enumerate(itertools.pairwise(cs)):
        for j, s in enumerate(part3.plan_edge(pt, z0, z1, part3.RHO_MAX_K * km)):
            s.update({"edge": edge0 + k, "idx": j})
            segs.append(s)
    return segs


def jsonl(l: int) -> Path:
    return RES3 / f"band_PB_l{l}_segments.jsonl"


def _init_retry(name):
    from . import part3
    part3._worker_init(name)
    part3.K_MAX = 240                      # the registered K5 retry (K_MAX 160 -> 240)


def run(ls=None, workers: int = 4, retry: bool = False):
    import multiprocessing as mpr

    from . import part3
    ls = L_SMALL if ls is None else ls
    RES3.mkdir(parents=True, exist_ok=True)
    pt = Point("PB", 256)
    rho_min = part3.RHO_MIN_K * float(pt.km.mid()) / (2 if retry else 1)
    tasks = []
    for l in ls:
        done = set()
        if jsonl(l).exists():
            for line in jsonl(l).read_text().splitlines():
                r = json.loads(line)
                if r.get("ok"):
                    done.add(r["key"])
        cs, _, _, _ = corners(pt, l)
        for s in plan_path(pt, cs):
            if part3.seg_key(s) not in done:
                tasks.append(("PB", l, part3._ser_seg(s), rho_min))
    tasks.sort(key=lambda t: -(t[1] * 1000 + abs(float(Fraction_(t[2]["c"][0]))) / float(pt.km.mid())))
    t_start = time.time()
    rec = {"cmd": sys.argv, "point": "PB", "ls": ls, "workers": workers, "retry": retry, "nice": os.nice(0),
           "host": socket.gethostname(), "pid": os.getpid(), "commit": commit_tag(),
           "start_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t_start)), "tasks": len(tasks)}
    lpath = RES3 / f"launch_{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime(t_start))}.json"
    lpath.write_text(json.dumps(rec, indent=1) + "\n")
    print(f"{len(tasks)} segments to run for l = {ls} (retry = {retry})", flush=True)
    nfail = 0
    ctx = mpr.get_context("spawn")
    init = _init_retry if retry else part3._worker_init
    with ctx.Pool(workers, initializer=init, initargs=("PB",)) as pool:
        for r in pool.imap_unordered(part3._task, tasks):
            if retry:
                r["retry"] = True
            with open(jsonl(r["l"]), "a") as fh:
                fh.write(json.dumps(r) + "\n")
            if not r["ok"]:
                nfail += 1
                print("FAILED", r["l"], r["key"], r["error"], flush=True)
    rec.update({"end_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "wall_s": round(time.time() - t_start, 1),
                "failed": nfail})
    lpath.write_text(json.dumps(rec, indent=1) + "\n")


def Fraction_(s):
    from fractions import Fraction
    return Fraction(s)


def launch(workers: int = 4, retry: bool = False):
    RES3.mkdir(parents=True, exist_ok=True)
    stamp = time.strftime("%Y%m%dT%H%M%S")
    log = RES3 / f"run_{stamp}.log"
    cmd = ["nice", "-n", "10", sys.executable, "-m", "p2.beta_part3", "run", str(workers)] + (["retry"] if retry else [])
    with log.open("w") as fh:
        p = subprocess.Popen(cmd, cwd=SRC, stdout=fh, stderr=subprocess.STDOUT, start_new_session=True)
    rec = {"pid": p.pid, "log": log.name, "cmd": " ".join(cmd), "commit": commit_tag(), "host": platform.node()}
    (RES3 / f"launch_cmd_{stamp}.json").write_text(json.dumps(rec, indent=1) + "\n")
    return rec


def _launch_commits(l: int) -> list:
    out = set()
    for f in sorted(RES3.glob("launch_2*.json")):
        rec = json.loads(f.read_text())
        if l in rec.get("ls", []):
            out.add(rec.get("commit"))
    return sorted(out)


def assemble(l: int, write: bool = True) -> dict:
    from . import part3
    pt = Point("PB", 256)
    with precision(part3.PREC):
        cs, a, bp, b0 = corners(pt, l)
        segs = plan_path(pt, cs)
        res = {}
        if jsonl(l).exists():
            for line in jsonl(l).read_text().splitlines():
                r = json.loads(line)
                if r["ok"]:
                    res[r["key"]] = r["results"]
        missing = [part3.seg_key(s) for s in segs if part3.seg_key(s) not in res]
        if missing:
            return {"l": l, "complete": False, "missing": len(missing)}
        pieces = []
        for s in segs:
            sub = res[part3.seg_key(s)]
            cur = [str(v) for v in s["p0"]]
            for r in sub:
                if r["p0"] != cur:
                    raise ArithmeticError(f"coverage gap in {part3.seg_key(s)}")
                cur = r["p1"]
            if cur != [str(v) for v in s["p1"]]:
                raise ArithmeticError(f"coverage gap at the end of {part3.seg_key(s)}")
            pieces += sub
        total = sum((arb(r["increment"]) for r in pieces), arb(0))
        wn = total / arb.pi()
        k = round(float(wn.mid()))
        cert = bool(abs(wn - k) < arb("0.25"))
        expected = 0 if l == 0 else -1
        e = arb(EDGE.numerator) / EDGE.denominator * pt.km
        out = {"claim": f"S6 part (iii) at PB, l = {l}: winding of W_l around the band box B_l "
                        "(zeros minus poles, with multiplicity)",
               "point": "PB", "l": l,
               "region": {"a_prime": str(a), "b_prime": str(bp), "b0": str(b0),
                          "a_source": "certificates/beta/part2/A_fixed_PB.json",
                          "b_prime_ge_e_kappa": bool(arb(bp) >= e), "b0_le_half_kappa": bool(arb(b0) <= pt.km / 2),
                          "omega0_in_interior": bool(arb(bp) > pt.kp and pt.kp > arb(b0))},
               "symmetry": "W(-conj w) = conj W(w): winding = (sum over the right half H)/pi",
               "sum_increments_over_H": total.str(30), "winding_ball": wn.str(20),
               "certified": cert, "winding": k if cert else None, "expected_winding": expected,
               "as_expected": bool(cert and k == expected),
               "accounting": ("Z = winding + P, P <= 1 (only omega0, order <= 1): l = 0 with omega_* inside forces "
                              "Z = P = 1; l >= 1 with winding -1 forces P = 1, Z = 0"),
               "segments": len(pieces), "K_max_used": max(r["K"] for r in pieces),
               "min_abs_W_on_H_lower": min(r["min_abs_enclosure_lower"] for r in pieces),
               "max_M_over_absWc": max(r["M_over_absWc"] for r in pieces),
               "cpu_s": round(sum(r["timing_s"]["cpu"] for r in pieces), 1),
               "method": {"prec_bits": part3.PREC, "theta_point": part3.THETA_POINT, "theta_box": part3.THETA_BOX,
                          "x": "1/2", "rho_max_kappa": part3.RHO_MAX_K, "target": part3.TARGET,
                          "retried_segments": sum(1 for r in pieces if r.get("split_from") is not None)},
               "computed_at": _launch_commits(l), "assembled_at": commit_tag(), "pieces": pieces}
    if write:
        write_json(CERT3 / f"band_PB_l{l}.json", out)
    return out


# ---------------------------------------------------------------------------------------------
# KATs and controls
# ---------------------------------------------------------------------------------------------
KAT_BOXES = [
    ("B2", 0, (-0.03, 0.03, -0.93, -0.89), 1, "omega_* only (top edge just above omega_*)"),
    ("B3", 0, (-0.03, 0.03, -0.86, -0.79), -1, "omega0 only (the pole)"),
    ("B3", 1, (-0.03, 0.03, -0.86, -0.79), -1, "omega0 only (the pole), l = 1"),
    ("B4", 0, (-0.03, 0.03, -0.93, -0.79), 0, "omega0 and omega_*"),
    ("B6", 0, (-0.05, 0.05, -1.90, -1.80), 1, "l = 0 NE overtone ~ -1.8552i kappa_- (below the edge; no pole inside)"),
    ("B6", 1, (6.9, 7.5, -1.80, -1.70), 1, "l = 1 PS mode ~ 7.2149 - 1.7533i kappa_- (below the edge)"),
]


def right_half_sum(pt: Point, l: int, cs: list, workers: int):
    from . import part3
    segs = plan_path(pt, cs)
    out = part3._pool_segments("PB", l, segs, workers, part3.RHO_MIN_K * float(pt.km.mid()))
    bad = [r for r in out if not r["ok"]]
    if bad:
        return None, bad[0]["error"], len(segs)
    with precision(part3.PREC):
        tot = sum((arb(x["increment"]) for r in out for x in r["results"]), arb(0))
    return tot, None, sum(len(r["results"]) for r in out)


def kat(workers: int = 4) -> dict:
    from . import part3
    pt = Point("PB", 256)
    out = {"commit": commit_tag(), "boxes": []}
    with Timer() as tm:
        for kid, l, bx, expect, what in KAT_BOXES:
            r = part3.rectangle_winding("PB", l, bx, workers=workers)
            r.update({"id": kid, "what": what, "expected": expect, "as_expected": bool(r.get("certified") and r.get("winding") == expect)})
            out["boxes"].append(r)
            print("KAT", kid, l, bx, r.get("winding"), r["as_expected"], flush=True)
        # NB1: the l = 0 band box with its bottom edge at -0.905 kappa_- (above omega_*): winding -1
        with precision(part3.PREC):
            bot = part3._dy_up(0.905 * float(pt.km.mid()))
            cs, a, _, b0 = corners(pt, 0, bottom=bot)
        tot, err, n = right_half_sum(pt, 0, cs, workers)
        with precision(part3.PREC):
            wn = (tot / arb.pi()) if tot is not None else None
            k = round(float(wn.mid())) if wn is not None else None
            certified = bool(wn is not None and abs(wn - k) < arb("0.25"))
        out["NB1"] = {"box": ["-a'_0", "a'_0", "-0.905 kappa_- (dyadic up)", "-b0"], "bottom": str(bot),
                      "winding_ball": wn.str(10) if wn is not None else None, "winding": k if certified else None,
                      "segments": n, "error": err, "expected": -1, "as_expected": bool(certified and k == -1)}
        print("NB1", out["NB1"]["winding"], out["NB1"]["as_expected"], flush=True)
        # NB2: an edge through omega_* (Im omega_* to 10 digits) must not certify
        r = part3.rectangle_winding("PB", 0, (-0.03, 0.03, -0.9093809790, -0.89), workers=workers)
        out["NB2"] = {"box": [-0.03, 0.03, -0.9093809790, -0.89], "certified": bool(r.get("certified")),
                      "error": r.get("error"), "fails_as_required": not r.get("certified")}
        print("NB2", out["NB2"]["fails_as_required"], (r.get("error") or "")[:80], flush=True)
        # left halves at l = 0 and 2 (the symmetry): the left-half sum equals the right-half sum
        out["left_halves"] = []
        for l in (0, 2):
            with precision(part3.PREC):
                csr, a, bp, b0 = corners(pt, l)
                csl = [(fmpq(0), -b0), (-a, -b0), (-a, -bp), (fmpq(0), -bp)]
            totl, err, n = right_half_sum(pt, l, csl, workers)
            st = json.loads((CERT3 / f"band_PB_l{l}.json").read_text()) if (CERT3 / f"band_PB_l{l}.json").exists() else None
            with precision(part3.PREC):
                ok = bool(totl is not None and st is not None and totl.overlaps(arb(st["sum_increments_over_H"])))
            out["left_halves"].append({"l": l, "sum_left": totl.str(20) if totl is not None else None,
                                       "sum_right_stored": st["sum_increments_over_H"] if st else None,
                                       "segments": n, "error": err, "overlap": ok})
            print("left half", l, ok, flush=True)
        # B1: replay of P2-K6's Krawczyk ball (krawczyk_PB_l0.json is not regenerated; A4)
        stored = json.loads((CERT.parent / "part3" / "krawczyk_PB_l0.json").read_text())
        kr = part3.krawczyk_z("PB", 0, ("0", "-0.9093809789948054555249"))
        with precision(256):
            ov = bool(arb(kr["enclosure_over_kappa_minus"]["im"]).overlaps(arb(stored["enclosure_over_kappa_minus"]["im"]))
                      and arb(kr["enclosure_over_kappa_minus"]["re"]).overlaps(arb(stored["enclosure_over_kappa_minus"]["re"])))
        out["B1"] = {"stored": stored["enclosure_over_kappa_minus"], "replayed": kr["enclosure_over_kappa_minus"],
                     "replayed_certified": kr["certified_unique_zero"], "contraction": kr["contraction"],
                     "overlap": ov, "pass": bool(kr["certified_unique_zero"] and ov)}
        print("B1", out["B1"]["pass"], flush=True)
    out["timing_s"] = {"wall": round(tm.wall, 1), "cpu": round(tm.cpu, 1)}
    out["pass"] = bool(all(b["as_expected"] for b in out["boxes"]) and out["NB1"]["as_expected"]
                       and out["NB2"]["fails_as_required"] and all(x["overlap"] for x in out["left_halves"])
                       and out["B1"]["pass"])
    write_json(CERT3 / "kat_controls_PB.json", out)
    return out


def replay(workers: int = 4, frac: float = 0.05, min_sample: int = 3) -> dict:
    """5% replay: a deterministic sample of every band certificate's pieces, recomputed with the stored K."""
    import multiprocessing as mpr

    from . import part3
    tasks = []
    per = {}
    for l in L_SMALL:
        st = json.loads((CERT3 / f"band_PB_l{l}.json").read_text())
        pieces = st["pieces"]
        n = len(pieces)
        k = max(min_sample, math.ceil(frac * n))
        idx = sorted({(j * 7919 + l * 104729) % n for j in range(k)})
        tasks += [("PB", l, pieces[j]) for j in idx]
        per[l] = len(idx)
    t0 = time.time()
    ctx = mpr.get_context("spawn")
    with ctx.Pool(workers, initializer=part3._worker_init, initargs=("PB",)) as pool:
        results = list(pool.imap(part3._replay_task, tasks))
    out = {"check": "5% replay of the band segments (stored K)", "sampled_per_l": per, "recomputed": len(results),
           "all_overlap": all(r["overlap"] for r in results), "results": results,
           "timing_s": {"wall": round(time.time() - t0, 1), "cpu_segments": round(sum(r["cpu_s"] for r in results), 1)},
           "commit": commit_tag()}
    out["pass"] = out["all_overlap"]
    write_json(CERT3 / "replay_PB.json", out)
    return out


def main(argv):
    cmd = argv[1] if len(argv) > 1 else "help"
    if cmd == "run":
        run(workers=int(argv[2]) if len(argv) > 2 else 4, retry=(len(argv) > 3 and argv[3] == "retry"))
    elif cmd == "launch":
        print(json.dumps(launch(int(argv[2]) if len(argv) > 2 else 4, len(argv) > 3 and argv[3] == "retry"), indent=1))
    elif cmd == "assemble":
        for l in L_SMALL:
            c = assemble(l)
            print(l, c.get("winding"), c.get("expected_winding"), c.get("as_expected"), c.get("segments"),
                  c.get("missing"), flush=True)
    elif cmd == "kat":
        c = kat(int(argv[2]) if len(argv) > 2 else 4)
        print("KAT/controls pass:", c["pass"])
    elif cmd == "replay":
        c = replay(int(argv[2]) if len(argv) > 2 else 4)
        print("replay pass:", c["pass"])
    else:
        print(__doc__)


if __name__ == "__main__":
    main(sys.argv)
