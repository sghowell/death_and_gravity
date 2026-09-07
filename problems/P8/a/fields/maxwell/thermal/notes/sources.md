# Primary sources, convention map and non-imports

Read directly on 2026-09-06. The actual-state specialization is supported
by primary field/stress results, while the continuous history bounds and
branch enclosures are independently derived in this child.

- [Fröb–Glavan–Meda, 2506.23193](https://arxiv.org/pdf/2506.23193),
  sections 3.1–3.2, propositions 3.3–3.4 and equations 3.21,3.24,3.25:
  the conformally thermal photon stress and exact background equation.
  Remark 3.2 explicitly cancels the derivative trace term by a finite
  curvature prescription. Our corresponding choice is the already
  available A.16 beta_M=0 member with no independent nonzero R^2
  gravitational coefficient; we do not erase that term for generic beta.
  Their coupling obeys kappa_source^2/2=kappa_here and their proper
  physical density is T_etaeta/a^2. Under the A.16 source-to-FK
  convention map their physical rho and pressure agree with our (2).
- The same paper's equation 3.2 contains scalar conformal factors in a
  gauge-potential notation. We do not copy that potential kernel as a
  coordinate field-strength covariance or settle its representation
  conventions here. The new proof directly constructs the two physical
  transverse occupations and uses A.16's gauge-invariant conformal map.
  This avoids importing an unnecessary ambiguous weight/index formula.
  The b_T used here is a length with explicit hbar, not the paper's
  dimensionally natural hbar=1 shorthand for inverse temperature.
- Corollary 3.5 and equations 3.27–3.30 of that paper discuss local
  solutions and a radiation approximation. We use neither its
  unspecified short interval nor that approximation to prove (11):
  the exact branch/domain is selected explicitly and the entire actual
  A.18 history is bounded by rational inequalities. Algebraic positive
  a^4 does not justify selecting an arbitrary root or crossing D=0.

The pinned [A.18 formulation](../../focusing/cosmology/FORMULATION.md),
[proof](../../focusing/cosmology/notes/proof.md) and
[closure assessment](../../focusing/cosmology/notes/closure.md) are relative
links from this notes directory. Its report, at
`../focusing/cosmology/certificates/cosmological-calibration.json` relative
to this child root, has SHA256
`b03767fb43c49a0f6c5ff34307355c4cd37a75930f61e07f4826c06101e4d23d`.
It and its full source-hashed A.17/A.16 lineage are replayed read-only.
A.16 supplies the physical anomaly normalization and source convention;
the independent engine reads its coefficients only after that transitive
pin is checked. No older scalar field/state or gamma value is reused.

The construction is an actual selected thermal state, not a theorem that
all Hadamard states solve the same geometry. It introduces no additional
classical radiation source. Its positive EED and finite endpoint are a
special exact-branch result, not a claimed new QEI-driven singularity or
all-beta existence theorem. Fundamental EFT validity at the endpoint,
interacting QED and observational fitting are not asserted. This
optional realization does not reopen the scoped objective already
achieved by the photon incompleteness/calibration theorem.
