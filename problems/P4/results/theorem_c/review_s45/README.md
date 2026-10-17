# S4-5 review checks (2026-10-17)

The scripts and outputs of the independent review
`docs/review-2026-10-17-p4-s4-5-adversarial.md` (§9 lists what each one checks).
They ran in a `git archive 75869b91` copy of `problems/P4`, with the gitignored
tube caches copied in:

    PYTHONPATH=<copy>/problems/P4/src python <script> [args]

- `cover_independent.py <copy>/problems/P4` → `cover_independent.txt`
- `s42_join.py <tube cache json>` → `s42_join.txt`
- `g1_centre.py` → `g1_centre.txt`
- every other script writes the JSON named after it, or the path given as its
  first argument

Captured stdout uses `.txt`, because the repository ignores `*.out` and
`*.log`.

`logs/` holds the stdout of the two end-to-end recomputes (sonic/s29 and
deep/deep-128). Both were bit-identical to the committed parts.
`merge_after_scratch_recompute.json` is the merge that followed: the strict rule
rejected both recomputed parts because they were written outside a git checkout
(review §8).
