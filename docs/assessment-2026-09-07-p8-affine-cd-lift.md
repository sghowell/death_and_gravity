# P8: an exact original-frame auxiliary lift, with UV work still open

This follows the [two scoped reduction-route verdicts](assessment-2026-09-07-p8-matching-route-verdicts.md).
Original P8 remains open. The completed scoped photon objective and
linear classification are unchanged; neither receives a new prerequisite.

## What has changed

[S6.37](../problems/P8/s6/matching/affine/FORMULATION.md) gives a different
kind of parent: the metric and an unrestricted affine connection are
independent. This is an exact classical auxiliary representation, not
a two-metric derivative expansion. Eliminating the connection reproduces
the complete original CD/M1 action on an open neighborhood of its entire
closed clock tube, with the physical metric and free chi unchanged.

The match includes all five quadratic Hessian operators, the curvature
coefficient, the original scalar F away from X=1, and the lower-order
boundary terms. It is not just a rolling-background fit or a match of
the highest-derivative equations. The nonzero coefficients that the
ordinary retained-order two-metric action could not supply are present.
The old singular Horndeski-map transformation is not used.

The construction uses the independent-connection action of
[Aoki and Shimada, version 2](https://arxiv.org/pdf/1806.02589v2), but
derives its connection equation and stationary action independently.
Writing x=-X_repo, h=(1+u^2)^3 and w=h-1-x, its central inverse is
p=sqrt(w/h)/2, with p^2=f/2 and source denominator Delta=1.
The remaining scalar coefficients are fixed by a regular linear x-ODE,
not by an integral over physical time. Its solution is analytic on the
open tube and has a uniform explicit bound.

## Why the full connection audit mattered

The calculation retains all 64 connection components. Removing exactly
four projective gauge directions leaves a 60-dimensional algebraic
Hessian, whose determinant is -2^52 p^72(8p^2-1)^3. Thus Delta=1 alone
does not establish invertibility: p^2=1/8 has three additional unsourced
zero modes even though the particular scalar forcing remains compatible.

This exceptional locus is outside the original closed tube. There
p^2 is in [9/40,11/40], 8p^2-1>=4/5, and the exact block inverses give
row-sum norm at most 880109/36000<25 in the specified scalar rest frame
and trace gauge. This is a whole-domain rational bound, not a grid.

The literal elimination also exposes a missing X in the paper's printed
lower-order Q2 bracket. Both the complete action calculation and a
separate pure-Palatini control require p-3X*p_X, not p-3p_X. The printed
version is retained as a negative control. An along-clock check would
miss the defect, so the open-domain action and scalar boundary matter.

All 64 original connection Euler equations vanish on the selected
solution. The functional chain rule then identifies the reduced metric,
clock and matter equations with the target's equations. The known CD/M1
solution therefore lifts with its original physical geometry and matter
backreaction. No new nonlinear-stability claim is inferred.

## What this does not settle

The connection is auxiliary in this action. An indefinite algebraic
Hessian is not a propagating ghost diagnosis, and its inverse bound is
not a heavy mass or cutoff. Exact classical rewriting alone does not
pass the complete adopted B gate, much less the vacuum and finite-gravity
V/G gates. The scoped result must not be relabeled a UV completion.

The next discriminating task is a specified kinetic extension which
actually acts on the unresolved connection quotient, followed by its
full constraint, spectrum and quantitative reduction analysis. Merely
adding a healthy spectator field would leave the target's UV problem
unchanged. A different vacuum extension must also be named and checked;
the original coefficient tube is not silently changed to obtain one.

## Verification

The frozen report records 20 source hashes, 114 named exact identities
comprising 1617 scalar entries, 30 independent/interface proof checks
and 46 rejected-input controls. The ordinary checkpoint replay passes
208 tests in 122.18 seconds; the independent read-only certificate CLI
also passes. Neither uses the broad-regression GCD adapter.

The complete P8 regression, with no test-directory exclusions, passes
3803 tests in 688.08 seconds. Its opt-in exact Gaussian-GCD runner first
checks 128 normalized reference tuples and records 6509 exact descents
and 4846 ordinary-domain fallbacks. Both random seeds are zero, and the
diagnosed host faulthandler timer plugin is disabled. A separate
adversarial action/domain review found no unresolved defect; these are
symbolic tests and written proofs, not a Lean formalization or external
peer review.
