"""S3 review: the author's gitignored per-segment working files against the committed certificates.

Reads (read only) <author worktree>/problems/P2/results/part3/PB_l{l}_segments.jsonl and checks:
one ok record per planned segment, no failures, no splits, no duplicates, and each record equal to
the corresponding piece of certificates/part3/winding_PB_l{l}.json.
Usage: python compare_author_jsonl.py [author_results_dir]
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT = HERE.parent / "part3"
DEFAULT = Path("/Users/seanhowell/dev/research/death_and_gravity/.claude/worktrees/p2-s3/problems/P2/results/part3")


def main():
    A = Path(sys.argv[1]) if len(sys.argv) > 1 else DEFAULT
    total, issues = 0, []
    for l in range(20):
        lines = [json.loads(x) for x in (A / f"PB_l{l}_segments.jsonl").read_text().splitlines()]
        keys = [r["key"] for r in lines]
        if len(keys) != len(set(keys)):
            issues.append(f"l={l}: duplicate keys")
        if not all(r["ok"] for r in lines):
            issues.append(f"l={l}: failed records")
        if not all(len(r["results"]) == 1 for r in lines):
            issues.append(f"l={l}: split segments")
        res = {r["key"]: r["results"][0] for r in lines}
        cert = json.loads((CERT / f"winding_PB_l{l}.json").read_text())
        for p in cert["pieces"]:
            if res.get(p["key"]) != p:
                issues.append(f"l={l}: piece {p['key']} differs")
        if len(res) != len(cert["pieces"]):
            issues.append(f"l={l}: {len(res)} records vs {len(cert['pieces'])} pieces")
        total += len(lines)
    print(json.dumps({"records": total, "issues": issues, "ok": not issues}))


if __name__ == "__main__":
    main()
