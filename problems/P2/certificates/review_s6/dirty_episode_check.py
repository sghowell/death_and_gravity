"""S6 review, check A2: the +dirty episode (read-only comparison).

Reads, without writing, the author's git-ignored L run files in the author's worktree
(.claude/worktrees/p2-s6/problems/P2/results/beta/part1/L.w0.first-pass-dirty.jsonl and L.w0.jsonl)
and compares the 178 first-pass (+dirty) rows with the committed certificate
certificates/beta/part1/L_PB.json.  -> dirty_episode_check_out.json
"""
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
AUTH = HERE.parents[4] / "p2-s6" / "problems" / "P2" / "results" / "beta" / "part1"
st = json.loads((HERE.parent / "beta" / "part1" / "L_PB.json").read_text())
committed = {tuple(x["box"]): x for x in st["leaves_detail"]}
dirty = [json.loads(x) for x in (AUTH / "L.w0.first-pass-dirty.jsonl").read_text().splitlines() if x.strip()]
final = [json.loads(x) for x in (AUTH / "L.w0.jsonl").read_text().splitlines() if x.strip()]
same = 0
diff = []
for r in dirty:
    for lf in r["leaves"]:
        c = committed.get(tuple(lf["box"]))
        if c == lf:
            same += 1
        else:
            diff.append(r["grid_box"])
fin_commits = sorted({r["commit"] for r in final})
fin_keys = [tuple(r["grid_box"]) for r in final]
out = {"dirty_rows": len(dirty), "dirty_commit_tags": sorted({r["commit"] for r in dirty}),
       "dirty_leaves_identical_to_committed": same, "dirty_leaves_different": diff[:10],
       "final_rows": len(final), "final_unique_grid_boxes": len(set(fin_keys)), "final_commit_tags": fin_commits,
       "dirty_grid_boxes_recomputed_in_final": sum(1 for r in dirty if tuple(r["grid_box"]) in set(fin_keys)),
       "read_only": "no file in the author's worktree was written"}
(HERE / "dirty_episode_check_out.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps(out, indent=1))
