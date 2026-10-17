#!/usr/bin/env python3
"""S4-6 review: the numbers notes/s4-6-centre.md quotes, against the data.

  * section 5.2's per-piece table (bound to the printed digits, boxes, evals,
    cpu seconds) against parts/centre_t14 and centre_t224;
  * section 5.1's CPU sum;
  * section 5.3's merge claims against results/theorem_c/tc_commutator.json;
  * section 5.4's sampling claims against results/theorem_c/centre_sampling.json;
  * that s4-5-assembly.md's S4-6 edits replaced lines in place (git diff
    hunks with equal old/new line counts, except the appended section 9).

usage (repository root): python problems/P4/results/theorem_c/review_s46/notes_check.py [out.json]
"""
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", "..", "..", "..", ".."))
P4 = os.path.join(ROOT, "problems", "P4")
TC = os.path.join(P4, "results", "theorem_c")


def tol(txt):
    """Half a unit in the last printed place of a decimal like -0.0126 or -5.91e-4."""
    m, _, e = txt.lower().partition("e")
    dec = len(m.split(".")[1]) if "." in m else 0
    return 0.5 * 10.0 ** (int(e or 0) - dec) * (1 + 1e-9)


def main(out):
    rec = dict(schema="review-s46-notes-check/1")
    note = open(os.path.join(P4, "notes", "s4-6-centre.md")).read()
    rows = re.findall(r"^\| (c\d\d) \| \[([^]]+)\] \| (\S+) \| (\d+) \| (\d+) \| (\d+) \| (\S+) \| "
                      r"(\d+) \| (\d+) \|$", note, re.M)
    mism = []
    cpu14 = cpu224 = build = 0.0
    for r in rows:
        pid = r[0]
        d14 = json.load(open(os.path.join(TC, "parts/centre_t14/centre/%s.json" % pid)))
        d224 = json.load(open(os.path.join(TC, "parts/centre_t224/centre/%s.json" % pid)))
        b14, b224 = float(r[2]), float(r[6])
        checks = dict(
            bound14=abs(d14["certified_bound"] - b14) <= tol(r[2]),
            boxes14=d14["result"]["boxes"] == int(r[3]),
            evals14=d14["result"]["evals"] == int(r[4]),
            cpu14=round(d14["timing"]["certify_cpu_s"]) == int(r[5]),
            bound224=abs(d224["certified_bound"] - b224) <= tol(r[6]),
            boxes224=d224["result"]["boxes"] == int(r[7]),
            cpu224=round(d224["timing"]["certify_cpu_s"]) == int(r[8]))
        if not all(checks.values()):
            mism.append(dict(piece=pid, row=r, part14=[d14["certified_bound"],
                                                       d14["result"]["boxes"],
                                                       d14["result"]["evals"],
                                                       d14["timing"]["certify_cpu_s"]],
                             part224=[d224["certified_bound"], d224["result"]["boxes"],
                                      d224["timing"]["certify_cpu_s"]],
                             failed=[k for k, v in checks.items() if not v]))
        cpu14 += d14["timing"]["certify_cpu_s"]
        cpu224 += d224["timing"]["certify_cpu_s"]
        build += d14["timing"]["build_cpu_s"]
    rec["table"] = dict(rows=len(rows), mismatches=mism)
    rec["cpu"] = dict(build=build, t224=cpu224, t14=cpu14,
                      core_h=(build + cpu224 + cpu14) / 3600)
    m = json.load(open(os.path.join(TC, "tc_commutator.json")))
    cs = m["centre_substitution"]
    worst = {}
    for e in m["detail"]:
        reg = e["region"]
        b = e.get("certified_bound")
        if b is None and e.get("layer"):
            continue
        if b is not None and (reg not in worst or b > worst[reg][0]):
            worst[reg] = (b, e["id"])
    rec["merge"] = dict(
        status=m["status"], ok=m["ok"], certified=m["strip_exclusion_certified"],
        substituted_layer=m["substitution"]["count"], substituted_centre=cs["count"],
        missing=len(m["missing"]), failed=len(m["failed"]),
        manual=len(m["substitution"]["manual_review"]) + len(cs["manual_review"]),
        rejected=len(m["substitution"]["rejected"]) + len(cs["rejected"]),
        blocked=len(m["substitution"]["blocked"]) + len(cs["blocked"]),
        pieces=m["cover"]["pieces"], cover_errors=m["cover"]["errors"],
        gaps=[m["scope"]["tau_le_tstar_uncovered_x"], m["scope"]["tau_gt_tstar_uncovered_x"]],
        open_at_centre_end=m["scope"]["open_at_centre_end"],
        checks_valid=m["checks_valid_for_code"], conflicts=m["code_conflicts"],
        accepted=m["accepted_code_variants"], identity_valid=cs["identity"]["valid"],
        second_mode=cs["second_mode"]["states"],
        worst_by_region={k: [v[0], v[1]] for k, v in sorted(worst.items())},
        git_head=m["provenance"]["git_head"])
    s = json.load(open(os.path.join(TC, "centre_sampling.json")))
    per = {k: v["max_direct_mid"] for k, v in s["pieces"].items()}
    rec["sampling"] = dict(points=s["points"], failures=s["failures"], ok=s["ok"],
                           sup_direct_mid=s["sup_direct_mid"],
                           max_abs_graded_minus_direct_mid=s["max_abs_graded_minus_direct_mid"],
                           direct_finite=sum(v["direct_finite"] for v in s["pieces"].values()),
                           graded_defined=sum(v["graded_defined"] for v in s["pieces"].values()),
                           tight=sum(v["direct_mid_tight"] for v in s["pieces"].values()),
                           max_graded_rad=max(v["max_graded_rad"] for v in s["pieces"].values()),
                           max_direct_mid_c27_c28_c29=[per["centre/c27"], per["centre/c28"],
                                                       per["centre/c29"]],
                           git_head=s["code"]["git_head"], dirty=s["code"]["git_dirty_src"])
    diff = subprocess.run(["git", "-C", ROOT, "diff", "main", "--",
                           "problems/P4/notes/s4-5-assembly.md"],
                          capture_output=True, text=True).stdout
    hunks = re.findall(r"^@@ -(\d+),(\d+) \+(\d+),(\d+) @@", diff, re.M)
    rec["assembly_hunks"] = [list(map(int, h)) for h in hunks]
    rec["assembly_lines_in_place"] = all(int(h[0]) == int(h[2]) and int(h[1]) == int(h[3])
                                         for h in hunks[:-1])
    json.dump(rec, open(out, "w"), indent=1)
    print(json.dumps(rec, indent=1)[:6000])


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "notes_check.json"))
