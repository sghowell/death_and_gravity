# Finite-mass first two-gauge cut

This immutable child of S6.128 bounds the entire momentum
Taylor remainder of the finite-mass fermion box, then sews
both factors of the first two-gauge cut with the correct
crossed boundary signs.

At the unchanged prospective mF=10^200 boundary, the relative
error at s=3/2 is below 5 times 10^-194. This retains a nonzero
perturbative cut inside the old massive scalar disc.
The result is not an all-orders amplitude or P8 closure.

See FORMULATION.md and notes/ for the proof and scope.
The read-only entry point is
python -m p8_vacuum_finite_mass_gauge_cut.verify --check,
with all P8 src directories on the existing import path.
No verifier writes its own expected report. Independent
science, every-field mutation tests and exact source hashes
separate finite algebra from the written analytic arguments.
