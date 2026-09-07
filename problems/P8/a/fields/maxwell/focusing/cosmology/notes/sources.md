# Source and hypothesis audit

The following primary sources and immutable repository inputs were
checked directly on 2026-09-06. The new arithmetic has a proof in this
child; it is not attributed to a paper that did not contain it.

- [Fewster–Pfenning, gr-qc/0303106](https://arxiv.org/pdf/gr-qc/0303106),
  theorem V.1 and section VI.A, equations 97–102: the physical spin-one
  all-Hadamard positive-type inequality and the Maxwell polarization
  coefficient. The paper explicitly removes its topological restrictions
  for Minkowski spacetime. A.16 supplies the local conformal-strip
  argument, not a claim that every strip state extends to global
  Minkowski. The present child imports that established local result.
- [Herzog–Huang, 1301.5002](https://arxiv.org/pdf/1301.5002), equation 23:
  type-D-free conformal reference stress. The A.16 primary audit derives
  its FK signature/curvature map and adds a separately defined conserved
  finite R^2 variation. No derivative tensor is imported with an unchecked
  source convention, and beta_M is not replaced by the scalar gamma.
- [Fewster–Kontou, 1907.13604](https://arxiv.org/pdf/1907.13604), equations
  11–12 and section 4.1.2: the timelike index and past-extended sampler
  structure. A.17 retains an exact geometric past square and proves the
  global FLRW conclusion itself. Neither its geometric history nor the
  new conditional future bounds are consequences of this literature QEI.

Relative to this notes directory, the immediate theorem inputs are
[A.17 formulation](../../FORMULATION.md), [proof](../../notes/proof.md)
and [source audit](../../notes/sources.md). Its report is
`../certificates/history-focusing.json` relative to the A.18 child root,
with SHA256
`22cfba547368b420767bbe87068c5b45fb556992c9e37a389091525f3e2ef62a`.
The verifier checks that hash, replays its full source-hashed lineage,
and only then reads its cubic moments and the linked A.16 photon EED
coefficient polynomial. All old source and certificate files remain
unchanged. Reading A.16 through this transitive pin does not introduce
an unpinned alternative field coefficient.

The power-law family, its continuous p/x bounds, the positive C3 margins,
the finite-beta affine split and the physical source/clock equations are
derived here. The radiation/dust labels mean the classical geometric
comparison p=1/2 and p=2/3, not an actual quantum-state/SEE construction.
For p between them the corresponding classical perfect-fluid w is
2/(3p)-1; this follows by substituting the power law into the classical
zero-Lambda Friedmann and continuity equations. No such fluid is added
to the photon theorem unless its source and energy assumptions are
separately specified.

No observations, numerical Hubble measurements or choice of photon beta_M
are imported. This calibration proves sufficient cosmological-scale
constants in a stated theorem, not empirical validity of its hypotheses
or an unrestricted completion claim. The distinction from the literal
original objective is made in [the closure assessment](closure.md).
