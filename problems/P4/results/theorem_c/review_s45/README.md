# S4-5 review checks (2026-10-17)

The scripts and outputs of the independent review
`docs/review-2026-10-17-p4-s4-5-adversarial.md` (§9 lists what each one checks).
They ran in a `git archive 75869b91` copy of `problems/P4`, with the gitignored
tube caches copied in:

    PYTHONPATH=<copy>/problems/P4/src python <script> [args]

- `cover_independent.py <copy>/problems/P4` → `cover_independent.out`
- `s42_join.py <tube cache json>` → `s42_join.out`
- `g1_centre.py` → `g1_centre.out`
- every other script writes the JSON named after it, or the path given as its
  first argument

`logs/` holds the two end-to-end recomputes (sonic/s29 and deep/deep-128). Both
were bit-identical to the committed parts, and the strict rule rejected both
because they were written outside a git checkout.

`merge_after_scratch_recompute.json` from the same run is not kept; review §8
summarizes it.
