"""S3 review: full rigorous replay of the author's pipeline (p2.part3.cauchy_segment, unchanged) for
every stored segment of the chosen l, with K re-derived (K=None), in this review worktree.

Compares each recomputed increment with the stored ball (overlap), checks that K and the piece
count are reproduced, and re-sums the winding from the recomputed increments.  Resumable: one
JSONL line per finished segment in logs/replay_l{l}.jsonl.  Writes replay_author.json.
Usage (cwd anywhere):  nice -n 10 python replay_author.py 1,2,19 --workers 3
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mpr
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parents[1] / "src"
sys.path.insert(0, str(SRC))

from flint import arb  # noqa: E402

from p2 import part3  # noqa: E402
from p2.arbutil import git_commit, precision  # noqa: E402
from p2.params import Point  # noqa: E402

CERT = HERE.parent / "part3"
LOGS = HERE / "logs"
_PT = {}


def _init():
    _PT["pt"] = Point("PB", 256)


def _job(args):
    l, r = args
    pt = _PT["pt"]
    seg = {"edge": r["edge"], "idx": r["idx"], "p0": [part3._fq(x) for x in r["p0"]],
           "p1": [part3._fq(x) for x in r["p1"]], "c": [part3._fq(x) for x in r["c"]],
           "rho": part3._fq(r["rho"])}
    t0 = time.process_time()
    with precision(part3.PREC):
        try:
            rr = part3.cauchy_segment(pt, l, seg, K=None)
        except ArithmeticError as ex:
            return {"l": l, "key": r["key"], "ok": False, "error": str(ex)}
        ov = bool(arb(rr["increment"]).overlaps(arb(r["increment"])))
    return {"l": l, "key": r["key"], "ok": True, "increment": rr["increment"], "stored": r["increment"],
            "overlap": ov, "K": rr["K"], "K_stored": r["K"], "pieces": rr["pieces"], "pieces_stored": r["pieces"],
            "M_upper": rr["M_upper"], "M_stored": r["M_upper"], "cpu": round(time.process_time() - t0, 2)}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ls")
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    ls = [int(x) for x in a.ls.split(",")]
    LOGS.mkdir(exist_ok=True)
    tasks = []
    stored = {}
    for l in ls:
        st = json.loads((CERT / f"winding_PB_l{l}.json").read_text())
        stored[l] = st
        f = LOGS / f"replay_l{l}.jsonl"
        done = set()
        if f.exists():
            for line in f.read_text().splitlines():
                rec = json.loads(line)
                if rec.get("ok"):
                    done.add(rec["key"])
        tasks += [(l, r) for r in st["pieces"] if r["key"] not in done]
    tasks.sort(key=lambda t: -(t[0] * 1000 + float(part3._fq(t[1]["c"][0]))))
    print(f"{len(tasks)} segments to replay; commit {git_commit()}", flush=True)
    ctx = mpr.get_context("spawn")
    with ctx.Pool(a.workers, initializer=_init) as pool:
        for rec in pool.imap_unordered(_job, tasks):
            with open(LOGS / f"replay_l{rec['l']}.jsonl", "a") as fh:
                fh.write(json.dumps(rec) + "\n")
            if not rec["ok"] or not rec["overlap"]:
                print("PROBLEM", rec, flush=True)
    # summary
    outp = HERE / "replay_author.json"
    out = json.loads(outp.read_text()) if outp.exists() else {"per_l": {}}
    out["code_commit"] = git_commit()
    for l in ls:
        recs = {}
        for line in (LOGS / f"replay_l{l}.jsonl").read_text().splitlines():
            rec = json.loads(line)
            recs[rec["key"]] = rec
        pieces = stored[l]["pieces"]
        missing = [r["key"] for r in pieces if r["key"] not in recs or not recs[r["key"]]["ok"]]
        with precision(part3.PREC):
            total = sum((arb(recs[r["key"]]["increment"]) for r in pieces if r["key"] in recs and recs[r["key"]]["ok"]), arb(0))
            wn = total / arb.pi()
            k = round(float(wn.mid()))
            cert = bool(not missing and abs(wn - k) < arb("0.25"))
        out["per_l"][str(l)] = {
            "segments": len(pieces), "replayed": len(pieces) - len(missing), "missing_or_failed": missing,
            "all_overlap": all(recs[r["key"]]["overlap"] for r in pieces if r["key"] in recs and recs[r["key"]]["ok"]),
            "K_reproduced": all(recs[r["key"]]["K"] == recs[r["key"]]["K_stored"] for r in pieces if r["key"] in recs and recs[r["key"]]["ok"]),
            "pieces_reproduced": all(recs[r["key"]]["pieces"] == recs[r["key"]]["pieces_stored"] for r in pieces if r["key"] in recs and recs[r["key"]]["ok"]),
            "winding_ball": wn.str(10), "certified": cert, "winding": k if cert else None,
            "stored_winding": stored[l]["winding"],
            "cpu_s": round(sum(recs[r["key"]].get("cpu", 0) for r in pieces if r["key"] in recs), 1)}
        print(l, out["per_l"][str(l)], flush=True)
    outp.write_text(json.dumps(out, indent=1) + "\n")


if __name__ == "__main__":
    main()
