# S5 — Perturbative control, starting with D-only M0

Opened 2026-09-04 following authorization to continue P8. This is an extension
of the published linear checkpoint `417724fc8a1086636dd960825ec36145b1e36999`,
not a relaxation or redefinition of its acceptance criteria. S5 is not yet
complete. No UV or strong-coupling verdict is inferred from S0--S4.

## Target and gates

Start with the exact D-only M0 witness. Restore a physical bounce duration
tau>0 and reduced Planck mass M>0. Determine whether a local fluctuation EFT
has an energy window above the background scales throughout the bounce
and both infinite tails. A failed test of this one witness is not a D-row
exclusion. M1 and UV dispersion/positivity are subsequent, separate gates.

1. **S5.1, nonlinear preparation and audit:** derive an invertible nonlinear
   field chart, verify its boundary terms and Legendre transform, recover
   the old linear constraint coefficients independently, and eliminate the
   lapse through fourth invariant order without gamma denominators. Restore
   physical scales and identify which remaining constraints prevent a cutoff
   conclusion. Compactify time in local background units and bound the
   normalized invariant coefficients in both tails. A scoped exact replay
   can certify these identities and bounds.
2. **S5.2, physical interaction reduction:** solve the spatial momentum
   constraints and retain generic momenta and scalar/tensor polarizations.
   Derive the canonical cubic and quartic interactions in covering regular
   charts. Track time-dependent canonical boundaries and operator-redefinition
   terms. A one-mode cubic spatial average is not an interaction calculation.
3. **S5.3, control test:** compute or conservatively bound the physical
   interactions, including cubic exchange and quartic contacts and mixed
   channels. State frequency, angular/IR and adiabatic assumptions. Compare
   to curvature and coefficient-variation scales, not |H| alone. Prove
   all-time coverage; report a hierarchy only for the sector/order actually
   checked. All-orders control needs a separate power-counting/convergence
   argument. Cancellations can invalidate single-vertex cutoff estimates.
4. **S5.4, later applicability:** assess boost breaking, a scattering domain,
   massless graviton poles and UV assumptions before using any positivity
   bound. No P9 positivity prior is imported. No UV completion assertion.

## Evidence and preservation

Self-review with independently derived algebra is not external peer review.
Written calculus/implicit-function lemmas are not Lean FORMALIZED results.
Use the existing local locked environment, exact arithmetic and regression
tests; no paid compute or external publication is needed for these steps.

New work lives only under `s5/` so S0--S4's source-glob hashes still replay.
The S5 report pins its own inputs and the prior classification certificate.
If an actual old-proof defect is found, record and correct its claim status;
preserving a checkpoint never licenses preserving a known false claim.

The initial deliverable is a lapse-reduced Hamiltonian in exact metric and
momentum invariants, **not** a fully reduced scattering action. Finite lapse
coefficients at the bounce are necessary preparation, not weak-coupling
evidence. The lapse Hessian tends to zero in the tails, which must be treated
with physical normalization rather than a compact-interval argument.
The implemented background-unit compactification supplies uniform algebraic
lapse solvability and coefficient bounds, not a canonical interaction bound.
