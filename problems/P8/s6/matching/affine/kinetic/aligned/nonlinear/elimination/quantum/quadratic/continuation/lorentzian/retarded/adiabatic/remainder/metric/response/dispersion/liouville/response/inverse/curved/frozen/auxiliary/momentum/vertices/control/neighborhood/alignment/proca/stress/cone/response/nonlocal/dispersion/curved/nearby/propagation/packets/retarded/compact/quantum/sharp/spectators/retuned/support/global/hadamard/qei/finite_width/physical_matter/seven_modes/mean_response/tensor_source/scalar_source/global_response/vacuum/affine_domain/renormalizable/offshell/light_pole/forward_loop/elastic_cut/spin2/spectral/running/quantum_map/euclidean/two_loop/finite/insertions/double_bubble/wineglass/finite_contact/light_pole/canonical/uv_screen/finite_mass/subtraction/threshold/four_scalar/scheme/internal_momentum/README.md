# Global internal-momentum control for the selected one-loop kernel

P8 S6.134 proves an unbounded complex-half-plane bound for the fully
on-shell-subtracted fermion insertion, retains a global spacelike
logarithmic envelope, and checks all real momenta in one named bubble
routing. Both one-loop sectors give a fractional Euclidean inverse defect
below 10^-202 on a specified finite reference-energy window.

The window is not a justified cutoff. This is an ingredient for subsequent
loop bounds, not a complete two-loop calculation or original P8 closure.

See FORMULATION.md and all six notes for the analytic arguments and limits.
The read-only verifier recursively rebuilds its frozen S6.133 ancestor,
checks every source hash and all report fields, and does not save reports.

Run the repository runtime wrapper from the repository root:

    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_global_one_loop_insertions.verify

It grants interpreter recursion and trusted exact-integer formatting
allowances without altering the scientific SymPy operations. The ordinary
checkpoint tests use the same unmodified library. The repository full
regression's separately documented exact adapter is not native evidence.
