"""P2-3 (part iii, S3): validated zero counts of W_l on the compact boxes R_l.

Region.  R'_l = [-a_l, a_l] x [-b, e] with exact dyadic corners (absolute units):
  * a_l is the certified P2-2 threshold (certificates/lg/lg_{X}_l{l}.json, field W_l_exact);
  * b = beta = (1/2 + delta) kappa_- rounded UP to 2^-30, so R'_l contains {|Re w| <= a_l, -beta <= Im w <= e};
  * e = eta = kappa_-/10 rounded DOWN to 2^-30 (0 < e <= eta; Im w > 0 is covered by P2-0).
W_l is analytic on a neighbourhood of R'_l: its poles lie on E_+ u E_c = {-i n kappa_+, -i n kappa_c :
n >= 1}, all with Im w <= -kappa_+ < -b.

Symmetry.  W_l(-conj w) = conj W_l(w) (the phase-factored equations of zsol have coefficients that
are conjugated by w -> -conj w, and Z_h(r_h) = 1), so the winding number of W_l around the boundary
of R'_l equals (Delta_H arg W_l)/pi, where H is the right half of the boundary, counterclockwise:
the bottom edge from (0, -b) to (a_l, -b), the right edge up to (a_l, e), the top edge back to (0, e).

Discs (Cauchy interpolation, S1 winding.py, with the sup bound from zsol).  H is covered by
consecutive segments [p_j, p_{j+1}] with exact dyadic endpoints; each lies within rho_j/2 of its
centre c_j.  On the disc |w - c_j| <= rho_j (inside the box [c_j +- rho_j]^2, which is checked to
avoid E_+ u E_c):
  * M_j >= sup |W_l| on the box (zsol.wronskian_z box evaluation);
  * K point values on |w - c_j| = rho_j/2 give the coefficients W_k with the aliasing error
    M x^K/((1 - x^K) rho^k), x = 1/2;
  * W_l(c + d) lies in P(d) + disc(M (|d|/rho)^K/(1 - |d|/rho)) for |d| <= rho/2.
Each segment is cut into pieces until every piece's enclosure excludes 0; the increment over a
piece is the principal Arg of the ratio of its end enclosures (certified |increment| < pi).  The
sum over H divided by pi must lie within 1/4 of an integer: that integer is the winding number.
"""
from __future__ import annotations

import argparse
import itertools
import json
import math
import os
import socket
import sys
import time
from fractions import Fraction
from pathlib import Path

from flint import acb, arb, fmpq

from .arbutil import git_commit, precision
from .params import Point
from . import zsol

PREC = 128
THETA_POINT = 3.0
THETA_BOX = 0.25
K_MIN, K_MAX = 24, 160
TARGET = 1e-6           # truncation + aliasing <= TARGET * min |W_c| (design; the enclosures decide)
RHO_MAX_K = 1.0         # largest disc radius, in units of kappa_-
RHO_MIN_K = 1.0 / 256   # smallest disc radius before a segment is declared not certifiable
PIECE_K = 0.0625        # initial piece length, in units of kappa_-
MAX_SPLIT = 14
DEN = 2 ** 30           # dyadic grid for positions and radii

ROOT = Path(__file__).resolve().parents[2]
CERT_DIR = ROOT / "certificates" / "part3"
RES_DIR = ROOT / "results" / "part3"
LG_DIR = ROOT / "certificates" / "lg"


# ---------------------------------------------------------------------------------------------
# geometry

def _dy_down(x: float) -> fmpq:
    return fmpq(math.floor(x * DEN), DEN)


def _dy_up(x: float) -> fmpq:
    return fmpq(math.ceil(x * DEN), DEN)


def _fq(s: str) -> fmpq:
    if "/" in s:
        n, d = s.split("/")
        return fmpq(int(n), int(d))
    return fmpq(int(s))


def region(pt: Point, l: int) -> dict:
    """Exact dyadic a_l, b >= beta, 0 < e <= eta (absolute units)."""
    st = json.loads((LG_DIR / f"lg_{pt.name}_l{l}.json").read_text())
    a = _fq(st["W_l_exact"])
    b = _dy_up(float(pt.beta.upper()))
    e = _dy_down(float(pt.eta.lower()))
    if not (arb(b) >= pt.beta and arb(e) <= pt.eta and arb(e) > 0 and arb(b) < pt.kp):
        raise ArithmeticError("region corners not certified")
    return {"a": a, "b": b, "e": e, "a_source": f"certificates/lg/lg_{pt.name}_l{l}.json"}


def poles_ok(pt: Point, c: acb, rho: arb, nmax: int = 12) -> bool:
    """Certify that the box [c +- rho]^2 (corner distance sqrt(2) rho) avoids E_+ u E_c."""
    s = arb(2).sqrt() * rho
    for kap in (pt.kp, pt.kc):
        for n in range(1, nmax + 1):
            if not bool(abs(c - acb(0, -n * kap)) > s):
                return False
    # poles with n > nmax lie below -nmax kappa_+ and are farther than the n = nmax ones from any
    # c with Im c > -nmax kappa_+ (checked here)
    return bool(c.imag > -nmax * pt.kp)


def _pole_dist_f(pt: Point, x: float, y: float) -> float:
    d = 1e300
    for kap in (float(pt.kp.mid()), float(pt.kc.mid())):
        for n in range(1, 13):
            d = min(d, math.hypot(x, y + n * kap))
    return d


def plan_edge(pt: Point, z0: tuple, z1: tuple, rho_max: float) -> list:
    """Segments covering the straight edge z0 -> z1 (exact (re, im) fmpq pairs; axis-parallel)."""
    horiz = z0[1] == z1[1]
    if horiz:
        s0, s1, fixed = z0[0], z1[0], z0[1]
    else:
        s0, s1, fixed = z0[1], z1[1], z0[0]
    sign = 1 if s1 > s0 else -1
    segs = []
    p = s0
    while (p - s1) * sign < 0:
        rho = rho_max
        while True:
            q = p + sign * _dy_down(rho)
            if (q - s1) * sign > 0:
                q = s1
            cpos = (p + q) / 2
            cx, cy = (float(cpos), float(fixed)) if horiz else (float(fixed), float(cpos))
            half = abs(float(q - p)) / 2
            rr = max(rho, 2 * half)
            if _pole_dist_f(pt, cx, cy) > 1.5 * math.sqrt(2) * rr:
                break
            rho /= 2
        rho_q = _dy_up(max(rho, 2 * half))
        if horiz:
            segs.append({"p0": [p, fixed], "p1": [q, fixed], "c": [cpos, fixed], "rho": rho_q})
        else:
            segs.append({"p0": [fixed, p], "p1": [fixed, q], "c": [fixed, cpos], "rho": rho_q})
        p = q
    return segs


def plan(pt: Point, l: int, rho_max_k: float = RHO_MAX_K) -> list:
    """The covering of H (right half of the boundary of R'_l), as a list of segment dicts."""
    reg = region(pt, l)
    a, b, e = reg["a"], reg["b"], reg["e"]
    km = float(pt.km.mid())
    rho_max = rho_max_k * km
    corners = [(fmpq(0), -b), (a, -b), (a, e), (fmpq(0), e)]
    segs = []
    for k, (z0, z1) in enumerate(itertools.pairwise(corners)):
        for j, s in enumerate(plan_edge(pt, z0, z1, rho_max)):
            s.update({"edge": k, "idx": j})
            segs.append(s)
    return segs


def plan_rectangle(pt: Point, x0: fmpq, x1: fmpq, y0: fmpq, y1: fmpq, rho_max_k: float = RHO_MAX_K) -> list:
    """The covering of the full boundary of [x0, x1] x [y0, y1] (counterclockwise; controls)."""
    km = float(pt.km.mid())
    corners = [(x0, y0), (x1, y0), (x1, y1), (x0, y1), (x0, y0)]
    segs = []
    for k, (z0, z1) in enumerate(itertools.pairwise(corners)):
        for j, s in enumerate(plan_edge(pt, z0, z1, rho_max_k * km)):
            s.update({"edge": k, "idx": j})
            segs.append(s)
    return segs


def seg_key(s: dict) -> str:
    return f"{s['edge']}:{s['p0'][0]}:{s['p0'][1]}:{s['p1'][0]}:{s['p1'][1]}"


# ---------------------------------------------------------------------------------------------
# one segment

def _acb_q(z) -> acb:
    return acb(arb(z[0]), arb(z[1]))


def _poly_eval(coef, d):
    acc = acb(0)
    for c in reversed(coef):
        acc = acc * d + c
    return acc


def _seg_ball(z0: acb, z1: acb) -> acb:
    mid = (z0 + z1) / 2
    hr = abs((z1 - z0).real) / 2
    hi = abs((z1 - z0).imag) / 2
    return acb(mid.real + arb(0, 1) * hr, mid.imag + arb(0, 1) * hi)


def _excludes_zero(z: acb) -> bool:
    return not z.real.contains(0) or not z.imag.contains(0)


def cauchy_segment(pt: Point, l: int, seg: dict, K: int | None = None) -> dict:
    """Arg increment of W_l along one segment (all rigorous; raises ArithmeticError on failure)."""
    t0, c0 = time.time(), time.process_time()
    c = _acb_q(seg["c"])
    rho = arb(seg["rho"])
    p0, p1 = _acb_q(seg["p0"]), _acb_q(seg["p1"])
    if not poles_ok(pt, c, rho):
        raise ArithmeticError("box meets E_+ or E_c")
    for z in (p0, p1):
        if not bool(abs(z - c) <= rho / 2):
            raise ArithmeticError("segment not within rho/2 of the centre")
    box = acb(c.real + arb(0, 1) * rho, c.imag + arb(0, 1) * rho)
    _, Wb = zsol.wronskian_z(pt, l, c, box, theta=THETA_BOX)
    Wc = zsol.wronskian_z(pt, l, c, theta=THETA_POINT)[0]
    M = arb(abs(Wb).upper())
    mlow = abs(Wc).lower()
    if not (mlow > 0) or not M.is_finite():
        raise ArithmeticError("centre value not separated from 0, or M not finite")
    ratio = float((M / arb(mlow)).upper())
    if K is None:
        need = math.log2(4 * ratio / TARGET)
        K = max(K_MIN, 8 * math.ceil(need / 8))
    if K > K_MAX:
        raise ArithmeticError(f"K = {K} > K_MAX (M/|W_c| = {ratio:.2e})")
    x = arb(1) / 2
    rp = x * rho
    vals = []
    for j in range(K):
        th = 2 * arb.pi() * j / K
        w = c + rp * acb(th.cos(), th.sin())
        vals.append(zsol.wronskian_z(pt, l, w, theta=THETA_POINT)[0])
    alias = M * x ** K / (1 - x ** K)
    coef = []
    for k in range(K):
        Sk = acb(0)
        for j in range(K):
            th = -2 * arb.pi() * ((j * k) % K) / K
            Sk += vals[j] * acb(th.cos(), th.sin())
        Sk = Sk / K
        err = alias / rho ** k
        coef.append(Sk / rp ** k + acb(arb(0, 1) * err, arb(0, 1) * err))

    def encl(z):
        d = z - c
        s_ = arb(abs(d).upper())
        if not (s_ < rho):
            raise ArithmeticError("evaluation point outside the disc")
        rem = M * (s_ / rho) ** K / (1 - s_ / rho)
        return _poly_eval(coef, d) + acb(arb(0, 1) * rem, arb(0, 1) * rem)

    length = float(abs(p1 - p0).upper())
    n_init = max(1, math.ceil(length / (PIECE_K * float(pt.km.mid()))))
    total = arb(0)
    pieces = 0
    min_abs = None
    max_rel = 0.0
    stack = [(p0 + (p1 - p0) * k / n_init, p0 + (p1 - p0) * (k + 1) / n_init, 0) for k in range(n_init)][::-1]
    while stack:
        z0, z1, depth = stack.pop()
        E = encl(_seg_ball(z0, z1))
        if not _excludes_zero(E):
            if depth >= MAX_SPLIT:
                raise ArithmeticError("piece enclosure contains 0 at the maximal depth")
            zm = (z0 + z1) / 2
            stack.append((zm, z1, depth + 1))
            stack.append((z0, zm, depth + 1))
            continue
        inc = (encl(z1) / encl(z0)).arg()
        if not (abs(inc) < arb.pi()):
            raise ArithmeticError("piece increment not certified below pi")
        total += inc
        pieces += 1
        al = abs(E).lower()
        min_abs = al if min_abs is None else min(min_abs, al)
        max_rel = max(max_rel, float(max(E.real.rad(), E.imag.rad()).upper()) / float(abs(E).mid()))
    return {"key": seg_key(seg), "edge": seg["edge"], "idx": seg["idx"],
            "p0": [str(v) for v in seg["p0"]], "p1": [str(v) for v in seg["p1"]],
            "c": [str(v) for v in seg["c"]], "rho": str(seg["rho"]),
            "increment": total.str(40), "K": K, "M_upper": float(M.upper()),
            "absWc_lower": float(mlow), "M_over_absWc": ratio, "pieces": pieces,
            "min_abs_enclosure_lower": float(min_abs) if min_abs is not None else None,
            "max_enclosure_rel_radius": max_rel,
            "timing_s": {"wall": round(time.time() - t0, 2), "cpu": round(time.process_time() - c0, 2)}}


def segment_with_split(pt: Point, l: int, seg: dict, rho_min: float) -> list:
    """cauchy_segment, halving the segment (and its radius) on failure, down to rho_min."""
    try:
        return [cauchy_segment(pt, l, seg)]
    except ArithmeticError as ex:
        if float(seg["rho"]) / 2 < rho_min:
            raise ArithmeticError(f"segment {seg_key(seg)} not certifiable: {ex}")
        p0, p1 = seg["p0"], seg["p1"]
        pm = [(p0[0] + p1[0]) / 2, (p0[1] + p1[1]) / 2]
        out = []
        for a_, b_ in ((p0, pm), (pm, p1)):
            child = dict(seg)
            child.update({"p0": a_, "p1": b_, "c": [(a_[0] + b_[0]) / 2, (a_[1] + b_[1]) / 2],
                          "rho": seg["rho"] / 2, "split_from": seg_key(seg), "split_reason": str(ex)[:120]})
            out += segment_with_split(pt, l, child, rho_min)
        return out


# ---------------------------------------------------------------------------------------------
# runs (resumable; one JSONL line per finished segment)

_PT = {}


def _worker_init(name):
    _PT["pt"] = Point(name, 256)


def _task(args):
    name, l, seg_s, rho_min = args
    seg = _de_seg(seg_s)
    pt = _PT.get("pt") or Point(name, 256)
    with precision(PREC):
        try:
            res = segment_with_split(pt, l, seg, rho_min)
            return {"ok": True, "l": l, "key": seg_key(seg), "results": res}
        except ArithmeticError as ex:
            return {"ok": False, "l": l, "key": seg_key(seg), "error": str(ex)}


def _ser_seg(s: dict) -> dict:
    return {k: ([str(x) for x in v] if isinstance(v, list) else (str(v) if isinstance(v, fmpq) else v))
            for k, v in s.items()}


def _de_seg(s: dict) -> dict:
    out = dict(s)
    for k in ("p0", "p1", "c"):
        out[k] = [_fq(x) for x in s[k]]
    out["rho"] = _fq(s["rho"])
    return out


_RES = {"dir": RES_DIR}


def jsonl_path(name: str, l: int) -> Path:
    return _RES["dir"] / f"{name}_l{l}_segments.jsonl"


def run(name: str, ls: list, workers: int = 4, rho_min_k: float = RHO_MIN_K, launch: dict | None = None):
    import multiprocessing as mpr
    _RES["dir"].mkdir(parents=True, exist_ok=True)
    pt = Point(name, 256)
    rho_min = rho_min_k * float(pt.km.mid())
    tasks = []
    for l in ls:
        done = set()
        f = jsonl_path(name, l)
        if f.exists():
            for line in f.read_text().splitlines():
                r = json.loads(line)
                if r.get("ok"):
                    done.add(r["key"])
        for s in plan(pt, l):
            if seg_key(s) not in done:
                tasks.append((name, l, _ser_seg(s), rho_min))
    # most expensive first (high l, high frequency)
    tasks.sort(key=lambda t: -(t[1] * 1000 + float(_fq(t[2]["c"][0])) / float(pt.km.mid())))
    print(f"{len(tasks)} segments to run for l = {ls}", flush=True)
    t_start = time.time()
    rec = {"cmd": sys.argv, "point": name, "ls": ls, "workers": workers, "nice": os.nice(0),
           "host": socket.gethostname(), "pid": os.getpid(), "commit": git_commit(),
           "start_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(t_start)),
           "tasks": len(tasks), "resdir": str(_RES["dir"])}
    if launch:
        rec.update(launch)
    lpath = _RES["dir"] / f"launch_{time.strftime('%Y%m%dT%H%M%SZ', time.gmtime(t_start))}.json"
    lpath.write_text(json.dumps(rec, indent=1) + "\n")
    nfail = 0
    ctx = mpr.get_context("spawn")
    with ctx.Pool(workers, initializer=_worker_init, initargs=(name,)) as pool:
        for r in pool.imap_unordered(_task, tasks):
            with open(jsonl_path(name, r["l"]), "a") as fh:
                fh.write(json.dumps(r) + "\n")
            if not r["ok"]:
                nfail += 1
                print("FAILED", r["l"], r["key"], r["error"], flush=True)
    rec.update({"end_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "wall_s": round(time.time() - t_start, 1),
                "failed": nfail})
    lpath.write_text(json.dumps(rec, indent=1) + "\n")


def _launch_commits(name: str, l: int) -> list:
    """Commits of the launch records (in the results directory) whose runs included l."""
    out = set()
    for f in sorted(_RES["dir"].glob("launch_*.json")):
        rec = json.loads(f.read_text())
        if rec.get("point") == name and l in rec.get("ls", []):
            out.add(rec.get("commit"))
    return sorted(out)


def assemble(name: str, l: int, write: bool = True) -> dict:
    """Check the coverage of H exactly, sum the increments, write the certificate."""
    pt = Point(name, 256)
    with precision(PREC):
        segs = plan(pt, l)
        res = {}
        for line in jsonl_path(name, l).read_text().splitlines():
            r = json.loads(line)
            if r["ok"]:
                res[r["key"]] = r["results"]
        missing = [seg_key(s) for s in segs if seg_key(s) not in res]
        if missing:
            return {"l": l, "complete": False, "missing": len(missing)}
        pieces = []
        for s in segs:
            sub = res[seg_key(s)]
            # the sub-segments must tile the planned segment in order, with exact shared endpoints
            cur = [str(v) for v in s["p0"]]
            for r in sub:
                if r["p0"] != cur:
                    raise ArithmeticError(f"coverage gap in {seg_key(s)}")
                cur = r["p1"]
            if cur != [str(v) for v in s["p1"]]:
                raise ArithmeticError(f"coverage gap at the end of {seg_key(s)}")
            pieces += sub
        total = sum((arb(r["increment"]) for r in pieces), arb(0))
        wn = total / arb.pi()
        k = round(float(wn.mid()))
        cert = bool(abs(wn - k) < arb("0.25"))
        reg = region(pt, l)
        out = {"claim": f"P2-3 at {name}, l = {l}: number of zeros of W_l in R'_l (with multiplicity)",
               "point": name, "l": l, "region": {"a": str(reg["a"]), "b": str(reg["b"]), "e": str(reg["e"]),
                                                 "a_source": reg["a_source"],
                                                 "b_ge_beta": True, "e_le_eta": True},
               "symmetry": "W(-conj w) = conj W(w): winding = (sum over the right half H)/pi",
               "sum_increments_over_H": total.str(30), "winding_ball": wn.str(20),
               "certified": cert, "winding": k if cert else None,
               "segments": len(pieces), "K_max_used": max(r["K"] for r in pieces),
               "min_abs_W_on_H_lower": min(r["min_abs_enclosure_lower"] for r in pieces),
               "max_M_over_absWc": max(r["M_over_absWc"] for r in pieces),
               "cpu_s": round(sum(r["timing_s"]["cpu"] for r in pieces), 1),
               "method": {"prec_bits": PREC, "theta_point": THETA_POINT, "theta_box": THETA_BOX,
                          "x": "1/2", "rho_max_kappa": RHO_MAX_K, "target": TARGET},
               "computed_at": _launch_commits(name, l), "assembled_at": git_commit(), "pieces": pieces}
    if write:
        CERT_DIR.mkdir(parents=True, exist_ok=True)
        (CERT_DIR / f"winding_{name}_l{l}.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


# ---------------------------------------------------------------------------------------------
# controls (KAT-style; notes/s3-preregistration.md section 4)

def _kq(pt: Point, x: float) -> fmpq:
    """Exact dyadic near x * kappa_- (positions of the control rectangles; any exact value is valid)."""
    return fmpq(round(x * float(pt.km.mid()) * DEN), DEN)


def _pool_segments(name: str, l: int, segs: list, workers: int, rho_min: float) -> list:
    import multiprocessing as mpr
    tasks = [(name, l, _ser_seg(s), rho_min) for s in segs]
    if workers <= 1:
        _worker_init(name)
        return [_task(t) for t in tasks]
    ctx = mpr.get_context("spawn")
    with ctx.Pool(workers, initializer=_worker_init, initargs=(name,)) as pool:
        return list(pool.imap(_task, tasks))


def rectangle_winding(name: str, l: int, box_k: tuple, workers: int = 1) -> dict:
    """Winding number of W_l around the full boundary of [x0, x1] x [y0, y1] (kappa_- units)."""
    pt = Point(name, 256)
    x0, x1, y0, y1 = (_kq(pt, v) for v in box_k)
    segs = plan_rectangle(pt, x0, x1, y0, y1)
    rho_min = RHO_MIN_K * float(pt.km.mid())
    out = _pool_segments(name, l, segs, workers, rho_min)
    bad = [r for r in out if not r["ok"]]
    if bad:
        return {"point": name, "l": l, "box_kappa": box_k, "certified": False, "error": bad[0]["error"]}
    with precision(PREC):
        pieces = [x for r in out for x in r["results"]]
        total = sum((arb(r["increment"]) for r in pieces), arb(0))
        wn = total / (2 * arb.pi())
        k = round(float(wn.mid()))
        cert = bool(abs(wn - k) < arb("0.25"))
    return {"point": name, "l": l, "box_kappa": list(box_k),
            "box_exact": [str(x0), str(x1), str(y0), str(y1)], "winding_ball": wn.str(15),
            "certified": cert, "winding": k if cert else None, "segments": len(pieces),
            "cpu_s": round(sum(r["timing_s"]["cpu"] for r in pieces), 1),
            "min_abs_W_lower": min(r["min_abs_enclosure_lower"] for r in pieces)}


def left_half_sum(name: str, l: int, workers: int = 1) -> dict:
    """Sum of the increments over the LEFT half of the boundary of R'_l (symmetry check)."""
    pt = Point(name, 256)
    reg = region(pt, l)
    a, b, e = reg["a"], reg["b"], reg["e"]
    km = float(pt.km.mid())
    corners = [(fmpq(0), e), (-a, e), (-a, -b), (fmpq(0), -b)]
    segs = []
    for k, (z0, z1) in enumerate(itertools.pairwise(corners)):
        for j, s_ in enumerate(plan_edge(pt, z0, z1, RHO_MAX_K * km)):
            s_.update({"edge": 10 + k, "idx": j})
            segs.append(s_)
    out = _pool_segments(name, l, segs, workers, RHO_MIN_K * km)
    bad = [r for r in out if not r["ok"]]
    if bad:
        return {"l": l, "certified": False, "error": bad[0]["error"]}
    with precision(PREC):
        total = sum((arb(x["increment"]) for r in out for x in r["results"]), arb(0))
    return {"l": l, "sum_left": total.str(30), "segments": sum(len(r["results"]) for r in out)}


def _arbw_reference(name: str, l: int, w_q: tuple, precs=(256, 512, 1024, 2048, 4096)):
    """S1's arbw.wronskian at the exact point w_q (fmpq pair), at the lowest precision giving a
    relative radius below 1e-10.  arbw's Frobenius stopping rule can end before its own tail
    condition N > 2|rho| + 1 holds at large |w| (it then raises); the minimum term count is raised
    to 2|rho| + 8 here, through a wrapper.  arbw itself is not modified."""
    from . import arbw
    orig = arbw.frobenius

    def frob(pt, which, w, L, t, min_terms=20):
        kap = pt.kp if which == "+" else pt.kc
        rho_abs = float((abs(acb(w)) / (2 * kap)).upper())
        return orig(pt, which, w, L, t, min_terms=max(min_terms, int(2 * rho_abs) + 8))

    arbw.frobenius = frob
    try:
        for prec in precs:
            pt = Point(name, prec)
            with precision(prec):
                w = acb(arb(w_q[0]), arb(w_q[1]))
                Wa = arbw.wronskian(pt, l, w)
                rr = float((abs(Wa).rad() / abs(Wa).mid()))
                if rr < 1e-10:
                    return Wa, prec, rr
        return Wa, precs[-1], rr
    finally:
        arbw.frobenius = orig


def zsol_vs_arbw(name: str, cases: list) -> list:
    """C0: the phase-factored W (zsol, 256 bits) against the S1 Frobenius/Taylor W (arbw) at exact
    points.  Pass: both balls have relative radius < 1e-10 and they overlap."""
    pt = Point(name, 256)
    out = []
    for l, re_k, im_k in cases:
        w_q = (_kq(pt, re_k), _kq(pt, im_k))
        with precision(256):
            w = acb(arb(w_q[0]), arb(w_q[1]))
            Wz = zsol.wronskian_z(pt, l, w, theta=THETA_POINT)[0]
            rz = float(abs(Wz).rad() / abs(Wz).mid())
        Wa, pa, ra = _arbw_reference(name, l, w_q)
        with precision(256):
            ov = bool(Wz.overlaps(Wa))
            rel = float((abs(acb(Wz.real.mid(), Wz.imag.mid()) - acb(Wa.real.mid(), Wa.imag.mid())) / abs(Wa).mid()).upper())
        out.append({"l": l, "w_kappa": [re_k, im_k], "w_exact": [str(w_q[0]), str(w_q[1])],
                    "overlap": ov, "relrad_zsol": rz, "relrad_arbw": ra, "arbw_prec": pa,
                    "meaningful": bool(rz < 1e-10 and ra < 1e-10), "pass": bool(ov and rz < 1e-10 and ra < 1e-10),
                    "W_zsol": Wz.str(15), "W_arbw": Wa.str(15), "rel_diff_mid": rel})
    return out


def d4_check(name: str = "PB") -> dict:
    """B1: W_0(0) contains 0 and the Cauchy coefficient W_1 at 0 contains i (r_+^2 + r_c^2)/(r_+ r_c)."""
    pt = Point(name, 256)
    with precision(PREC):
        W0 = zsol.wronskian_z(pt, 0, acb(0), theta=THETA_POINT)[0]
        rho = arb(_kq(pt, 0.25))
        c = acb(0)
        if not poles_ok(pt, c, rho):
            raise ArithmeticError("D4 disc meets a pole")
        box = acb(arb(0, 1) * rho, arb(0, 1) * rho)
        _, Wb = zsol.wronskian_z(pt, 0, c, box, theta=THETA_BOX)
        M = arb(abs(Wb).upper())
        K = 48
        x = arb(1) / 2
        rp = x * rho
        vals = [zsol.wronskian_z(pt, 0, rp * acb((2 * arb.pi() * j / K).cos(), (2 * arb.pi() * j / K).sin()),
                                 theta=THETA_POINT)[0] for j in range(K)]
        S1 = sum((vals[j] * acb((-2 * arb.pi() * j / K).cos(), (-2 * arb.pi() * j / K).sin()) for j in range(K)), acb(0)) / K
        err = M * x ** K / (1 - x ** K) / rho
        W1 = S1 / rp + acb(arb(0, 1) * err, arb(0, 1) * err)
        exact = acb(0, 1) * (pt.rp ** 2 + pt.rc ** 2) / (pt.rp * pt.rc)
    return {"W0_at_0": W0.str(10), "W0_contains_0": bool(W0.contains(0)), "W1": W1.str(20),
            "exact": exact.str(20), "W1_contains_exact": bool(W1.overlaps(exact)),
            "pass": bool(W0.contains(0) and W1.overlaps(exact))}


def krawczyk_z(name: str, l: int, z_k: tuple, rel_r: str = "1e-20", rel_R: str = "1e-4") -> dict:
    """The l-th Krawczyk unique-zero enclosure with the zsol W (as winding.krawczyk, S1 KAT-6)."""
    import mpmath as mp
    pt = Point(name, 256)
    with precision(256):
        km = pt.km
        mp.mp.dps = 80
        zc = acb(arb(z_k[0]), arb(z_k[1])) * acb(km)
        z0 = acb(zc.real.mid(), zc.imag.mid())
        absz = abs(z0)
        hfd = arb((absz * arb("1e-30")).mid())
        for _ in range(3):                                   # Newton on midpoints (design only)
            Wz = zsol.wronskian_z(pt, l, z0, theta=THETA_POINT)[0]
            d = (zsol.wronskian_z(pt, l, z0 + hfd, theta=THETA_POINT)[0]
                 - zsol.wronskian_z(pt, l, z0 - hfd, theta=THETA_POINT)[0]) / (2 * hfd)
            z0 = z0 - acb(Wz.real.mid(), Wz.imag.mid()) / acb(d.real.mid(), d.imag.mid())
            z0 = acb(z0.real.mid(), z0.imag.mid())
        r = arb((absz * arb(rel_r)).mid())
        R = arb((absz * arb(rel_R)).mid())
        s_corner = arb(2).sqrt() * r
        box_R = acb(z0.real + arb(0, 1) * R, z0.imag + arb(0, 1) * R)
        if not poles_ok(pt, z0, R):
            raise ArithmeticError("Krawczyk box meets a pole")
        _, WR = zsol.wronskian_z(pt, l, z0, box_R, theta=THETA_BOX)
        M = arb(abs(WR).upper())
        Wz0 = zsol.wronskian_z(pt, l, z0, theta=THETA_POINT)[0]
        fd = (zsol.wronskian_z(pt, l, z0 + hfd, theta=THETA_POINT)[0]
              - zsol.wronskian_z(pt, l, z0 - hfd, theta=THETA_POINT)[0]) / (2 * hfd)
        e_fd = M * hfd * hfd / (R ** 3 * (1 - (hfd / R) ** 2))
        if not bool(R > s_corner):
            raise ArithmeticError("Krawczyk: need R > sqrt(2) r")
        e_z = 2 * M * s_corner / (R - s_corner) ** 2
        dW = fd + acb(arb(0, 1), arb(0, 1)) * (e_fd + e_z)
        Y = 1 / acb(fd.real.mid(), fd.imag.mid())
        Zc = acb(arb(0, 1) * r, arb(0, 1) * r)
        Kz = z0 - Y * Wz0 + (1 - Y * dW) * Zc
        inside = bool(abs(Kz.real - z0.real) < r and abs(Kz.imag - z0.imag) < r)
        Kk = Kz / acb(km)
    return {"point": name, "l": l, "certified_unique_zero": inside,
            "enclosure_over_kappa_minus": {"re": Kk.real.str(30), "im": Kk.imag.str(30)},
            "Z_halfwidth_rel": rel_r, "R_rel": rel_R, "contraction": float(abs(1 - Y * dW).upper()),
            "M_upper": float(M.upper())}


def verify(name: str, ls: list, sample_frac: float = 0.05, min_sample: int = 3) -> list:
    """C7 replay: exact re-sum and coverage of every certificate, and a deterministic sample of
    segments recomputed (their increments must overlap the stored balls)."""
    probs = []
    pt = Point(name, 256)
    for l in ls:
        f = CERT_DIR / f"winding_{name}_l{l}.json"
        if not f.exists():
            probs.append(f"missing {f.name}")
            continue
        st = json.loads(f.read_text())
        with precision(PREC):
            pieces = st["pieces"]
            total = sum((arb(r["increment"]) for r in pieces), arb(0))
            if not total.overlaps(arb(st["sum_increments_over_H"])):
                probs.append(f"l={l}: re-sum differs")
            wn = total / arb.pi()
            if not (st["certified"] and abs(wn - st["winding"]) < arb("0.25")):
                probs.append(f"l={l}: winding not reproduced")
            # coverage of H, in order, with exact shared endpoints
            plan_segs = plan(pt, l)
            ends = [[str(v) for v in plan_segs[0]["p0"]]] + [[str(v) for v in s_["p1"]] for s_ in plan_segs]
            cur = ends[0]
            for r in pieces:
                if r["p0"] != cur:
                    probs.append(f"l={l}: coverage gap at {r['key']}")
                    break
                cur = r["p1"]
            if cur != ends[-1]:
                probs.append(f"l={l}: coverage does not end at (0, e)")
            n = len(pieces)
            k = max(min_sample, math.ceil(sample_frac * n))
            idx = sorted({(j * 7919 + l * 104729) % n for j in range(k)})
            for j in idx:
                r = pieces[j]
                seg = {"edge": r["edge"], "idx": r["idx"], "p0": [_fq(x) for x in r["p0"]],
                       "p1": [_fq(x) for x in r["p1"]], "c": [_fq(x) for x in r["c"]], "rho": _fq(r["rho"])}
                rr = cauchy_segment(pt, l, seg, K=r["K"])
                if not arb(rr["increment"]).overlaps(arb(r["increment"])):
                    probs.append(f"l={l}: segment {r['key']} increment not reproduced")
    return probs


C2_RECTS = [("PB", 1, (1.0, 12.0, -2.0, 0.1), 1, "P_B l=1 PS 7.2149-1.7533i"),
            ("PB", 2, (1.0, 17.0, -2.0, 0.1), 1, "P_B l=2 PS 12.1965-1.7079i"),
            ("PB", 0, (-0.05, 0.05, -0.95, -0.87), 1, "P_B l=0 NE -0.9094i"),
            ("PA", 1, (3.19, 3.29, -0.85, -0.75), 1, "P_A l=1 PS (CCDHJ18 Table I) 3.2426-0.7958i")]
C3_RECTS = [("PB", 1, (6.9, 7.5, -1.80, -1.743), 1, "top edge 0.0103 kappa_- above the l=1 PS mode"),
            ("PB", 1, (6.9, 7.5, -1.743, -1.70), 0, "bottom edge 0.0103 kappa_- above the l=1 PS mode")]


def c0_cases(name: str = "PB") -> list:
    pt = Point(name, 256)
    km = float(pt.km.mid())
    cases = []
    for l in (0, 1, 2, 5, 10, 19):
        a = float(region(pt, l)["a"]) / km
        cases += [(l, round(0.3 * a, 4), -0.55), (l, round(0.7 * a, 4), 0.1), (l, round(a, 4), -0.55), (l, 1.0, 0.1)]
    return cases


def symmetry_check(workers: int = 4) -> dict:
    """C4: the left half of the boundary of R'_2 computed explicitly against the certificate's right half."""
    lh = left_half_sum("PB", 2, workers)
    rf = CERT_DIR / "winding_PB_l2.json"
    right = json.loads(rf.read_text())["sum_increments_over_H"] if rf.exists() else None
    with precision(PREC):
        sym = bool(right is not None and "sum_left" in lh and arb(lh["sum_left"]).overlaps(arb(right)))
    out = {"left": lh, "right": right, "pass": sym, "commit": git_commit()}
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    (CERT_DIR / "symmetry_PB_l2.json").write_text(json.dumps(out, indent=1) + "\n")
    return out


def controls(workers: int = 4) -> dict:
    """C0-C3 and C6 of notes/s3-preregistration.md section 3 (the gate before production; C4, C5
    and C7 use the production certificates)."""
    t0, c0 = time.time(), time.process_time()
    out = {"commit": git_commit()}
    out["C0_zsol_vs_arbw"] = zsol_vs_arbw("PB", c0_cases("PB"))
    out["C0_pass"] = all(x["pass"] for x in out["C0_zsol_vs_arbw"])
    print("C0", out["C0_pass"], [(x["l"], x["arbw_prec"]) for x in out["C0_zsol_vs_arbw"]], flush=True)
    out["C1_D4"] = d4_check("PB")
    out["C2"] = []
    for name, l, box, exp, label in C2_RECTS:
        r = rectangle_winding(name, l, box, workers)
        r.update({"expected": exp, "label": label, "pass": bool(r.get("certified") and r.get("winding") == exp)})
        out["C2"].append(r)
        print("C2", label, r.get("winding"), r.get("winding_ball"), flush=True)
    out["C3"] = []
    for name, l, box, exp, label in C3_RECTS:
        r = rectangle_winding(name, l, box, workers)
        wrong = bool(r.get("certified") and r.get("winding") != exp)
        r.update({"expected": exp, "label": label, "wrong_count": wrong,
                  "pass": bool(r.get("certified") and r.get("winding") == exp)})
        out["C3"].append(r)
        print("C3", label, r.get("winding"), r.get("winding_ball"), flush=True)
    kr = krawczyk_z("PB", 0, ("0", "-0.9093809789948054555249"))
    k6 = json.loads((ROOT / "certificates" / "kat6_winding_krawczyk.json").read_text())["cases"][1]["krawczyk"]
    with precision(256):
        ov = bool(arb(kr["enclosure_over_kappa_minus"]["im"]).overlaps(arb(k6["enclosure_over_kappa_minus"]["im"])))
    kr["overlaps_KAT6"] = ov
    out["C6_krawczyk_l0"] = kr
    out["C6_pass"] = bool(kr["certified_unique_zero"] and ov)
    out["pass"] = bool(out["C0_pass"] and out["C1_D4"]["pass"] and all(x["pass"] for x in out["C2"])
                       and not any(x["wrong_count"] for x in out["C3"]) and out["C6_pass"])
    out["C3_all_certified"] = all(x["pass"] for x in out["C3"])
    out["timing_s"] = {"wall": round(time.time() - t0, 1), "cpu_parent": round(time.process_time() - c0, 1),
                       "cpu_segments": round(sum(x.get("cpu_s", 0) for x in out["C2"] + out["C3"]), 1)}
    CERT_DIR.mkdir(parents=True, exist_ok=True)
    (CERT_DIR / "controls_PB.json").write_text(json.dumps(out, indent=1) + "\n")
    kcert = dict(kr)
    kcert.update({"claim": "unique zero of W_0 at P_B in Z = z0 +- r (1 + i); the l = 0 near-extremal mode",
                  "KAT6_enclosure_im": k6["enclosure_over_kappa_minus"]["im"], "commit": out["commit"]})
    (CERT_DIR / "krawczyk_PB_l0.json").write_text(json.dumps(kcert, indent=1) + "\n")
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(prog="p2.part3")
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--point", default="PB")
    r.add_argument("--l", required=True, help="e.g. 0-19 or 0,1,2")
    r.add_argument("--workers", type=int, default=4)
    r.add_argument("--resdir", default=None)
    a = sub.add_parser("assemble")
    a.add_argument("--point", default="PB")
    a.add_argument("--l", required=True)
    a.add_argument("--resdir", default=None)
    a.add_argument("--no-write", action="store_true")
    v = sub.add_parser("verify")
    v.add_argument("--point", default="PB")
    v.add_argument("--l", required=True)
    v.add_argument("--sample", type=float, default=0.05)
    cc = sub.add_parser("controls")
    cc.add_argument("--workers", type=int, default=4)
    sy = sub.add_parser("symmetry")
    sy.add_argument("--workers", type=int, default=4)
    args = ap.parse_args(argv)
    if args.cmd == "symmetry":
        o = symmetry_check(args.workers)
        print("C4 symmetry pass:", o["pass"], flush=True)
        sys.exit(0 if o["pass"] else 1)
    if args.cmd == "controls":
        o = controls(args.workers)
        print("controls pass:", o["pass"], "C3 all certified:", o["C3_all_certified"], flush=True)
        sys.exit(0 if o["pass"] else 1)
    if getattr(args, "resdir", None):
        _RES["dir"] = Path(args.resdir)
    ls = []
    for part in (args.l.split(",") if getattr(args, "l", None) else []):
        if "-" in part:
            u, v = part.split("-")
            ls += list(range(int(u), int(v) + 1))
        else:
            ls.append(int(part))
    if args.cmd == "run":
        run(args.point, ls, args.workers)
    elif args.cmd == "verify":
        pr = verify(args.point, ls, args.sample)
        print("part3 replay:", "OK" if not pr else pr, flush=True)
        sys.exit(1 if pr else 0)
    else:
        for l in ls:
            o = assemble(args.point, l, write=not args.no_write)
            print(l, {k: o.get(k) for k in ("complete", "certified", "winding", "winding_ball", "segments", "cpu_s")}, flush=True)


if __name__ == "__main__":
    main()
