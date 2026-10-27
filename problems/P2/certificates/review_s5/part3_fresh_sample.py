#!/usr/bin/env python
"""S5 review: recompute a FRESH sample of S3 winding segments (P2-3B) with the author's code.

`python -m p2.part3 verify` always recomputes the same deterministic 5% sample (indices
(j*7919 + l*104729) mod n), which is the sample S3 itself committed (certificates/part3/replay_PB.json).
This script picks a different deterministic sample: 4 segments per l for l = 3..18 (the S3 review
replayed l = 0, 1, 2, 19 in full), avoiding the standard sample, and recomputes each with its stored K.
Each recomputed increment ball must overlap the stored one.  Nothing on disk is changed.

Run:  python part3_fresh_sample.py [workers]     (writes part3_fresh_sample_out.json)
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
P2 = HERE.parents[1]
sys.path.insert(0, str(P2 / "src"))

OUT = HERE / "part3_fresh_sample_out.json"


def pick(n: int, l: int, k: int = 4) -> list[int]:
    std = {(j * 7919 + l * 104729) % n for j in range(max(3, math.ceil(0.05 * n)))}
    out, j = [], 0
    while len(out) < k and j < 10 * n:
        i = (j * 104729 + l * 7919 + 31337) % n
        if i not in std and i not in out:
            out.append(i)
        j += 1
    return sorted(out)


def main(workers: int = 2):
    from p2 import part3
    t0 = time.time()
    tasks, chosen = [], {}
    for l in range(3, 19):
        st = json.loads((P2 / "certificates" / "part3" / f"winding_PB_l{l}.json").read_text())
        idx = pick(len(st["pieces"]), l)
        chosen[l] = idx
        tasks += [("PB", l, st["pieces"][i]) for i in idx]
    import multiprocessing as mpr
    ctx = mpr.get_context("spawn")
    with ctx.Pool(workers, initializer=part3._worker_init, initargs=("PB",)) as pool:
        res = list(pool.imap(part3._replay_task, tasks))
    out = {"tool": "part3_fresh_sample.py (S5 review; author's cauchy_segment, fresh sample)",
           "chosen_indices": chosen, "segments": len(res), "all_overlap": all(r["overlap"] for r in res),
           "results": res, "cpu_segments_s": round(sum(r["cpu_s"] for r in res), 1),
           "wall_s": round(time.time() - t0, 1)}
    OUT.write_text(json.dumps(out, indent=1) + "\n")
    print(f"fresh sample: {len(res)} segments, all overlap = {out['all_overlap']}, "
          f"cpu {out['cpu_segments_s']} s, wall {out['wall_s']} s")
    return 0 if out["all_overlap"] else 1


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 2))
