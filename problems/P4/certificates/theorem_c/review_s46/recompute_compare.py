#!/usr/bin/env python3
"""S4-6 review: the 60 centre parts recomputed in the review worktree against
the committed ones, and the merge on the recomputed parts.

Before running: `s46_centre.py work --force` has rewritten
results/theorem_c/parts/centre_t{14,224}/centre/*.json in the worktree.  This
script compares each with `git show HEAD:<path>` (every certificate field:
status, bounds, worst, boxes, evals, empty, depth histogram, the cover
parameters, ladder attempts, forms, engine, identity, inputs, scope, the
recorded source hashes and driver; only timing, env, dates and git head may
differ), then runs the real merge CLI
(`s43_parts.py merge --substitute-layer-t14 --substitute-centre`, which writes
results/theorem_c/tc_commutator.json), keeps a copy as merge_recomputed.json
and compares it with the committed tc_commutator.json (only the provenance's
git head and merge time may differ).  Restore the tracked files afterwards with
`git checkout -- problems/P4/results/theorem_c/parts/centre_t14
problems/P4/results/theorem_c/parts/centre_t224 problems/P4/results/theorem_c/tc_commutator.json`.

usage (repository root): python problems/P4/results/theorem_c/review_s46/recompute_compare.py
"""
import glob
import json
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
TC = os.path.join(ROOT, "problems", "P4", "results", "theorem_c")
PY = sys.executable

SAME = ("schema", "id", "region", "spec", "status", "tau_min", "scope", "certified_bound",
        "forms", "error", "engine", "inputs", "identity", "driver_sha256", "provenance")
RESULT_SAME = ("ok", "worst", "boxes", "evals", "empty", "fails", "depth_hist", "params",
               "certified_bound", "laurent", "norm", "scale", "relative_margin")


def head_json(rel):
    out = subprocess.run(["git", "-C", ROOT, "show", "HEAD:%s" % rel], capture_output=True,
                         check=True).stdout
    return json.loads(out)


def strip_attempts(a):
    return [(x["target"], x["ok"], x["evals"]) for x in a]


def main():
    rows, diffs = [], []
    for t in (14, 224):
        for p in sorted(glob.glob(os.path.join(TC, "parts", "centre_t%d" % t, "centre", "*.json"))):
            rel = os.path.relpath(p, ROOT)
            new, old = json.load(open(p)), head_json(rel)
            bad = [k for k in SAME if new.get(k) != old.get(k)]
            bad += ["result." + k for k in RESULT_SAME
                    if new["result"].get(k) != old["result"].get(k)]
            if strip_attempts(new["result"]["attempts"]) != strip_attempts(old["result"]["attempts"]):
                bad.append("result.attempts")
            if new["code"]["files"] != old["code"]["files"]:
                bad.append("code.files")
            prov = dict(git_head=new["code"]["git_head"], dirty=new["code"]["git_dirty_src"],
                        late=new["code"]["late_modules"], write_check=new["code"]["write_check"]["ok"])
            recomputed = new["timing"]["start"] != old["timing"]["start"]
            rows.append(dict(part=rel, recomputed=recomputed, identical=not bad, differing=bad,
                             bound=new["certified_bound"], boxes=new["result"]["boxes"],
                             evals=new["result"]["evals"],
                             cpu=new["timing"]["certify_cpu_s"] + (new["timing"]["build_cpu_s"]
                                                                   if t == 14 else 0.0),
                             provenance=prov))
            if bad or not recomputed:
                diffs.append(rows[-1])
    # the merge on the recomputed parts (writes the tracked tc_commutator.json)
    r = subprocess.run([PY, os.path.join(ROOT, "problems/P4/scripts/s43_parts.py"), "merge",
                        "--substitute-layer-t14", "--substitute-centre"],
                       cwd=ROOT, capture_output=True, text=True)
    out = os.path.join(TC, "tc_commutator.json")
    shutil.copy(out, os.path.join(HERE, "merge_recomputed.json"))
    new_m = json.load(open(out))
    old_m = head_json("problems/P4/results/theorem_c/tc_commutator.json")

    def walk(a, b, path=""):
        if isinstance(a, dict) and isinstance(b, dict):
            out_ = []
            for k in sorted(set(a) | set(b)):
                out_ += walk(a.get(k), b.get(k), path + "." + str(k))
            return out_
        if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
            out_ = []
            for i, (x, y) in enumerate(zip(a, b)):
                out_ += walk(x, y, path + "[%d]" % i)
            return out_
        return [] if a == b else [path]

    mdiff = walk(new_m, old_m)
    allowed = {".provenance.git_head", ".provenance.merged"}
    rec = dict(schema="review-s46-recompute-compare/1",
               parts=len(rows), recomputed=sum(1 for x in rows if x["recomputed"]),
               identical=sum(1 for x in rows if x["identical"]), differing=diffs,
               git_heads=sorted({x["provenance"]["git_head"] for x in rows}),
               all_clean=all(x["provenance"]["dirty"] is False and not x["provenance"]["late"]
                             and x["provenance"]["write_check"] for x in rows),
               cpu_core_h=sum(x["cpu"] for x in rows) / 3600,
               merge=dict(returncode=r.returncode, stdout=r.stdout.strip()[-600:],
                          status=new_m["status"], certified=new_m["strip_exclusion_certified"],
                          centre=new_m["centre_substitution"]["count"],
                          layer=new_m["substitution"]["count"],
                          diff_paths=mdiff, only_allowed=set(mdiff) <= allowed),
               rows=rows)
    rec["ok"] = (rec["parts"] == 60 and rec["recomputed"] == 60 and rec["identical"] == 60
                 and rec["all_clean"] and rec["merge"]["only_allowed"]
                 and rec["merge"]["certified"])
    json.dump(rec, open(os.path.join(HERE, "recompute_compare.json"), "w"), indent=1)
    print(json.dumps({k: v for k, v in rec.items() if k not in ("rows",)}, indent=1)[:3000])


if __name__ == "__main__":
    main()
