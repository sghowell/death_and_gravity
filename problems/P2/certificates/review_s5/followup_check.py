#!/usr/bin/env python
"""S5 review follow-up (2026-10-26): mechanical checks of the author's fix commit 3f2a17aa.

1. The P2-4 and P2-5 rows in notes/s5-report.md §6 are byte-identical to the review's §7 rows.
2. The .tex \\item[v1.3] in s5-report §7 equals the pre-review text (at 65af9e72) with exactly the review's
   two §7 corrections ((c), (d)) and the two approval/review dates filled in; the bibitem is unchanged.
3. The pre-registered plan (s5-plan.md at e2c75aa6) is a byte prefix of the current plan.
4. No residual scope wording ("not uniform", "up to the bifurcation sphere", "Reading of P2-5") in
   FORMULATION.md and the three S5 notes (verbatim HV17 quotes excepted).
5. 3f2a17aa changed nothing under src/, tests/, certificates/, docs/ or CLAIMS.md.
6. Writes followup_tex_item_v1.3.txt: the final \\item[v1.3], i.e. item 2's text with the one-word
   change of the review's §9 (5(i)): "one-ended Cauchy surface" -> "one-ended initial surface".

Run from anywhere inside the repository:  python followup_check.py   (writes followup_check_out.json)
"""
from __future__ import annotations

import hashlib
import json
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[3]
P2 = REPO / "problems" / "P2"
OUT = HERE / "followup_check_out.json"
FIX, REVIEWED_AT, PLAN_AT = "3f2a17aa", "65af9e72", "e2c75aa6"


def show(rev: str, rel: str) -> str:
    return subprocess.run(["git", "show", f"{rev}:{rel}"], cwd=REPO, capture_output=True, text=True,
                          check=True).stdout


def main():
    res = {}
    review = (REPO / "docs" / "review-2026-10-26-p2-s5-adversarial.md").read_text(encoding="utf-8")
    report = (P2 / "notes" / "s5-report.md").read_text(encoding="utf-8")

    # 1. rows
    sec6 = report.split("## 6.")[1].split("## 7.")[0].split("```")[1].split("\n")
    rows = {}
    for rid in ("| P2-4 |", "| P2-5 |"):
        a = [l for l in review.split("\n") if l.startswith(rid)]
        b = [l for l in sec6 if l.startswith(rid)]
        rows[rid.strip("| ")] = {"review_rows": len(a), "report_rows": len(b), "identical": a == b and len(a) == 1,
                                 "sha256": hashlib.sha256(b[0].encode()).hexdigest() if b else None}
    res["1_rows_identical"] = rows

    # 2. the .tex item
    def item(txt):
        return [l for l in txt.split("\n") if l.startswith("\\item[v1.3]")][0]

    def bib(txt):
        return [l for l in txt.split("\n") if l.startswith("\\bibitem{costafranzen2017}")]

    old_report = show(REVIEWED_AT, "problems/P2/notes/s5-report.md")
    old, new = item(old_report), item(report)
    reps = [
        ("lies in $H^{1.05-0}$, in particular in $H^1_{\\mathrm{loc}}$, at the Cauchy horizon.",
         "lies in $H^{1.05-0}$, in particular in $H^1_{\\mathrm{loc}}$, near the part of the Cauchy horizon in its "
         "domain of dependence; for smooth data on a complete Cauchy surface the same holds near every point of the "
         "Cauchy horizon other than its bifurcation sphere."),
        ("or uniform up to the bifurcation sphere of the Cauchy horizon", "or at the bifurcation sphere of the Cauchy horizon"),
        ("approved ⟨date⟩)", "approved October 26, 2026)"),
        ("independently reviewed ⟨date⟩;", "independently reviewed October 26, 2026;"),
    ]
    exp = old
    for a, b in reps:
        assert exp.count(a) == 1, a
        exp = exp.replace(a, b)
    res["2_tex_item"] = {"equals_old_plus_review_corrections_and_dates": new == exp,
                         "leading_placeholder_kept": new.startswith("\\item[v1.3] ⟨Month day⟩, 2026"),
                         "bibitem_unchanged": bib(old_report) == bib(report) and len(bib(report)) == 1}

    # 3. plan prefix
    plan_old = subprocess.run(["git", "show", f"{PLAN_AT}:problems/P2/notes/s5-plan.md"], cwd=REPO,
                              capture_output=True, check=True).stdout
    plan_new = (P2 / "notes" / "s5-plan.md").read_bytes()
    res["3_plan_prefix"] = {"pre_registered_text_is_prefix": plan_new.startswith(plan_old),
                            "bytes_registered": len(plan_old), "bytes_now": len(plan_new)}

    # 4. residual scope wording
    hits = []
    for rel in ("FORMULATION.md", "notes/s5-assembly.md", "notes/s5-report.md", "notes/s5-plan.md"):
        for i, line in enumerate((P2 / rel).read_text(encoding="utf-8").split("\n"), 1):
            low = line.lower()
            if ("not uniform" in low or "reading of p2-5" in low
                    or ("up to the bifurcation sphere" in low and "two-ended" not in low)):
                hits.append(f"{rel}:{i}")
    res["4_residual_scope_wording"] = {"hits": hits, "none": not hits}

    # 5. what the fix commit touched
    changed = subprocess.run(["git", "diff", "--name-only", REVIEWED_AT, FIX], cwd=REPO, capture_output=True,
                             text=True, check=True).stdout.split()
    bad = [c for c in changed if c.startswith(("problems/P2/src", "problems/P2/tests", "problems/P2/certificates",
                                                "docs/")) or c == "CLAIMS.md"]
    res["5_fix_commit_scope"] = {"changed": changed, "code_test_cert_docs_claims_untouched": not bad}

    # 6. the final .tex item with the one-word change (5(i))
    one = "on the one-ended Cauchy surface of \\cite{hintzvasy2017}"
    assert new.count(one) == 1
    final = new.replace(one, "on the one-ended initial surface of \\cite{hintzvasy2017}")
    (HERE / "followup_tex_item_v1.3.txt").write_text(final + "\n", encoding="utf-8")
    res["6_final_tex_item"] = {"file": "followup_tex_item_v1.3.txt", "change": "one word: Cauchy -> initial",
                               "sha256": hashlib.sha256(final.encode()).hexdigest()}

    ok = (all(v["identical"] for v in rows.values()) and res["2_tex_item"]["equals_old_plus_review_corrections_and_dates"]
          and res["2_tex_item"]["leading_placeholder_kept"] and res["2_tex_item"]["bibitem_unchanged"]
          and res["3_plan_prefix"]["pre_registered_text_is_prefix"] and res["4_residual_scope_wording"]["none"]
          and res["5_fix_commit_scope"]["code_test_cert_docs_claims_untouched"])
    res["pass"] = ok
    OUT.write_text(json.dumps(res, indent=1, ensure_ascii=False) + "\n", encoding="utf-8")
    for k, v in res.items():
        print(k, json.dumps(v, ensure_ascii=False)[:160])
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
