"""S6 review: replay the author's S6 certificates with the author's code, writing only into review_s6/.

Each job recomputes a committed certificate (or part of it) from the code at the review head and
compares the certified content with the committed file.  Nothing under certificates/beta or results/
is written: the jobs call the author's library functions directly (not the run/launch drivers).

    python replay_author.py epole | part2 | TA | controls | L | S+ | Sc | band SAMPLE_SEED NPIECES
(cwd: this directory; the script puts problems/P2/src on sys.path)  -> replay_<job>.json
"""
from __future__ import annotations

import hashlib
import json
import sys
import time
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT = HERE.parent / "beta"
sys.path.insert(0, str(HERE.parents[1] / "src"))          # problems/P2/src (the author's code under replay)


def _j(x):
    return json.loads(json.dumps(x, ensure_ascii=False))


def _strip(d, keys=("timing_s", "commit", "commits", "cpu_s_total", "timing", "seconds", "wall", "cpu")):
    if isinstance(d, dict):
        return {k: _strip(v, keys) for k, v in d.items() if k not in keys}
    if isinstance(d, list):
        return [_strip(v, keys) for v in d]
    return d


def _diff(a, b, path="", out=None, limit=40):
    out = [] if out is None else out
    if len(out) >= limit:
        return out
    if isinstance(a, dict) and isinstance(b, dict):
        for k in sorted(set(a) | set(b)):
            if k not in a or k not in b:
                out.append(f"{path}/{k}: present in only one")
            else:
                _diff(a[k], b[k], f"{path}/{k}", out, limit)
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            out.append(f"{path}: lengths {len(a)} vs {len(b)}")
        for i, (x, y) in enumerate(zip(a, b)):
            _diff(x, y, f"{path}[{i}]", out, limit)
    elif a != b:
        out.append(f"{path}: {str(a)[:80]} != {str(b)[:80]}")
    return out


def job_epole():
    from p2 import beta_epole
    new = _j(beta_epole.certificate())
    old = json.loads((CERT / "epole_PB.json").read_text())
    d = _diff(_strip(old), _strip(new))
    return {"identical_apart_from_timing_commit": not d, "differences": d, "new_commit": new.get("commit")}


def job_part2():
    from p2 import beta_part2
    new = _j(beta_part2.certificate())
    old = json.loads((CERT / "part2" / "A_fixed_PB.json").read_text())
    d = _diff(_strip(old), _strip(new))
    return {"identical_apart_from_timing_commit": not d, "differences": d, "pass": new["pass"],
            "a_prime": {r["l"]: r.get("a_prime") for r in new["rows"]}}


def job_TA():
    from p2 import beta_part1
    out = {}
    for name, fn in (("T", beta_part1.run_T), ("A", beta_part1.run_A)):
        new = _j(fn())
        old = json.loads((CERT / "part1" / f"{name}_PB.json").read_text())
        d = _diff(_strip(old), _strip(new))
        out[name] = {"identical_apart_from_timing_commit": not d, "differences": d, "pass": new["pass"]}
    new = _j(beta_part1.cover())
    old = json.loads((CERT / "part1" / "cover_PB.json").read_text())
    d = _diff(_strip(old), _strip(new))
    out["cover"] = {"identical_apart_from_commit": not d, "differences": d}
    return out


def job_controls():
    from p2 import beta_part1
    new = _j(beta_part1.controls())
    old = json.loads((CERT / "part1" / "controls_PB.json").read_text())
    d = _diff(_strip(old), _strip(new))
    return {"identical_apart_from_commit": not d, "differences": d, "pass": new["pass"]}


def job_L():
    from flint import arb  # noqa: F401
    from p2 import beta_common, beta_part1
    from p2 import lowfreq as LF
    from p2.arbutil import precision
    from p2.part1 import load_bulk
    W_bound, _ = load_bulk()
    st = json.loads((CERT / "part1" / "L_PB.json").read_text())
    stored = {tuple(x["box"]): x for x in st["leaves_detail"]}
    pt = beta_common.deep_point()
    same, diff, fails, n = 0, [], [], 0
    rows_sha = hashlib.sha256()
    t0 = time.time()
    with precision(LF.PREC):
        ch = beta_common.make_lchecker(pt, W_bound)
        for b in beta_part1.l_boxes():
            key = tuple(str(x) for x in b)
            leaves = _j(LF.certify_box_refined(ch, b, beta_common.L_MIN_BOX))
            n += 1
            for lf in leaves:
                rows_sha.update(json.dumps(lf, sort_keys=True).encode())
                if not lf["pass"]:
                    fails.append(key)
            if len(leaves) == 1 and tuple(leaves[0]["box"]) in stored and leaves[0] == stored[tuple(leaves[0]["box"])]:
                same += 1
            else:
                diff.append({"grid_box": key, "first": _diff(stored.get(key, {}), leaves[0] if leaves else {})[:5]})
    return {"grid_boxes": n, "identical_to_committed": same, "different": diff[:20], "n_different": len(diff),
            "failed": fails, "sha256_of_new_leaves": rows_sha.hexdigest(), "wall_s": round(time.time() - t0, 1)}


def job_S(side_name):
    from p2 import beta_common, beta_part1
    from p2 import throat as T
    from p2.arbutil import precision
    from p2.part1 import load_bulk
    beta_common.use_extended_s_grid()
    side = 1 if side_name == "+" else -1
    tag = "plus" if side > 0 else "c"
    W_bound, _ = load_bulk()
    st = json.loads((CERT / "part1" / f"S_PB_{tag}.json").read_text())
    stored = st["leaves"]
    pt = beta_common.deep_point()
    g = T.FGeo(pt)
    new = []
    t0 = time.time()
    with precision(T.PREC):
        for key in beta_part1.s_jobs():
            if key[0] == "closing":
                new += [T.certify_closing(pt, g, side, W_bound)]
            else:
                new += T.certify_zbox_refined(pt, g, side, Fraction(key[1]), Fraction(key[2]), W_bound)
    new = _j(new)
    exact = new == stored
    fm = []
    if not exact and len(new) == len(stored):
        for a, b in zip(stored, new):
            a2, b2 = dict(a), dict(b)
            fa, fb = a2.pop("float_margin", None), b2.pop("float_margin", None)
            if a2 != b2:
                fm.append({"leaf": a.get("zeta_over_OmegaP"), "diff": _diff(a2, b2)[:5]})
            elif fa != fb:
                fm.append({"leaf": a.get("zeta_over_OmegaP"), "float_margin_only": [fa, fb]})
    return {"leaves": len(new), "stored_leaves": len(stored), "exactly_identical": exact,
            "all_pass": all(x["pass"] for x in new), "non_identical": fm, "wall_s": round(time.time() - t0, 1)}


def job_band(seed: int, npieces: int):
    """recompute npieces band pieces chosen by a seeded hash different from the author's 5% sample,
    plus every piece within 0.2 kappa_- of omega0 (l = 0, 1) and every l = 0 piece within 0.05 kappa_- of omega_*"""
    import multiprocessing as mpr
    import random
    from p2 import part3
    from p2.params import Point
    pt = Point("PB", 256)
    km = float(pt.km.mid())
    kp = float(pt.kp.mid())
    rng = random.Random(seed)
    tasks = []
    author = json.loads((CERT / "part3" / "replay_PB.json").read_text())
    author_keys = {(r.get("l"), r.get("key")) for r in author["results"]}
    allp = []
    for l in range(20):
        pcs = json.loads((CERT / "part3" / f"band_PB_l{l}.json").read_text())["pieces"]
        for p in pcs:
            allp.append((l, p))
    special = []
    for l, p in allp:
        cx, cy = (float(Fraction(v)) for v in p["c"])
        if l in (0, 1) and abs(complex(cx, cy) - complex(0, -kp)) < 0.2 * km:
            special.append((l, p))
        elif l == 0 and abs(complex(cx, cy) - complex(0, -0.909381 * km)) < 0.05 * km:
            special.append((l, p))
    pool_ = [x for x in allp if (x[0], x[1]["key"]) not in author_keys and x not in special]
    pick = rng.sample(pool_, npieces)
    tasks = [("PB", l, p) for l, p in special + pick]
    t0 = time.time()
    ctx = mpr.get_context("spawn")
    with ctx.Pool(2, initializer=part3._worker_init, initargs=("PB",)) as pool:
        res = list(pool.imap(part3._replay_task, tasks))
    return {"seed": seed, "special_pieces_near_omega0_or_omega_star": len(special), "random_pieces": len(pick),
            "excluded_author_sample": len(author_keys), "recomputed": len(res),
            "all_overlap": all(r["overlap"] for r in res),
            "not_overlapping": [r for r in res if not r["overlap"]][:10],
            "per_l": {l: sum(1 for t in tasks if t[1] == l) for l in range(20)},
            "cpu_s_segments": round(sum(r.get("cpu_s", 0) for r in res), 1), "wall_s": round(time.time() - t0, 1),
            "results": [{"l": t[1], "key": t[2]["key"], "overlap": r["overlap"]} for t, r in zip(tasks, res)]}


def main(argv):
    job = argv[1]
    t0 = time.time()
    if job == "epole":
        out = job_epole()
    elif job == "part2":
        out = job_part2()
    elif job == "TA":
        out = job_TA()
    elif job == "controls":
        out = job_controls()
    elif job == "L":
        out = job_L()
    elif job in ("S+", "Sc"):
        out = job_S(job[1] if job == "S+" else "c")
    elif job == "band":
        out = job_band(int(argv[2]), int(argv[3]))
    else:
        raise SystemExit(__doc__)
    out["job"] = job
    out["total_wall_s"] = round(time.time() - t0, 1)
    name = {"S+": "S_plus", "Sc": "S_c"}.get(job, job)
    (HERE / f"replay_{name}.json").write_text(json.dumps(out, indent=1, ensure_ascii=False) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k not in ("results",)}, indent=1)[:3000])


if __name__ == "__main__":
    main(sys.argv)
