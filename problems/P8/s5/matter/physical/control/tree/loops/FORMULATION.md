# P8 S5.9.CD — Isolated M1 local one-loop audit

Adopted 2026-09-05. This new branch leaves the CD/M1 physical action,
S5.7 free-control and S5.8 finite-tree certificates unchanged. It pins and
replays S5.8 certificate SHA256
`839151b1ae103d4945b54ed703d9bb35e0aa4aa5ff3d173c20e7859767dd9787`.

## Claim and quantization scope

Integrate out only the one real, massless, minimally coupled free M1 scalar
chi in the original physical metric. The metric and DHOST clock remain
external classical fields. Determine the local one-loop ultraviolet pole,
its required counterterm class, and the size of its specified local
constant-scale logarithmic increment on CD. This is not the complete
coupled-theory beta function: mixed, clock and graviton loops, including
possible cancellations, are not calculated.

The strict quadratic-DHOST plus canonical-free-matter list is not closed
under this isolated determinant in its frozen field variables, even modulo
bulk total derivatives. The proof is off shell inside the permitted clock
tube, using a nonzero fourth-order metric Euler-Lagrange symbol. It does not
claim nonredundancy modulo arbitrary derivative field redefinitions or
equations of motion. A field redefinition may trade the new curvature
operators for other interactions; retaining the original physical matter
frame requires keeping the induced terms and observable map.

S6 already permits higher EFT operators. Accordingly this is a new matching
obligation, not an exclusion of the S6 admissibility class or a modification
of any old linear/tree verdict.

## Conventions and domains

P8 has signature +---, scalar curvature R=-6(Hdot+2H²), Einstein action
-M²R/2, and matter action +(partial chi)²/2. Keep hbar as the loop-counting
parameter, setting it to one only in the explicitly labelled numerical
example. Use the real-scalar Euclidean operator DE=-nablaE²+xi RE as a local
coefficient device; xi=0 is M1, and xi=1/6 is only a control, not a changed
witness. A2 is the s² coefficient after removing (4 pi s)^(-n/2).

The signed pole uses n=4-2 epsilon and GammaE=(hbar/2) Tr log(DE/mu²).
For the local continuation GammaE=-i GammaL|continued, dt=-i dtE and
gL|continued=-gE, curvature-squared coefficients satisfy cE=-cL. The table
in the proof fixes both pole and counterterm signs. Only the local
Lorentzian counterterm theorem is used; a global Wick-rotated CD geometry,
positive Euclidean determinant, infrared prescription or quantum state is
not assumed.

Curvature-squared decompositions and the zero bulk variations of Euler and
BoxR are four-dimensional statements with compactly supported variations.
Do not discard dimensionally continued Euler poles before renormalization:
evanescent finite/anomaly effects and boundary terms are outside this bound.
The selected BoxR representative comes directly from the general Laplace
coefficient and is not a claimed scheme-independent local observable.

All background estimates hold at every finite real u=t/tau, tau>0, with
a=(1+u²)², ell=tau sqrt(1+u²), and x=u/sqrt(1+u²). Tensor components are
lower orthonormal components. Bounds extend to the compact endpoint limits
|x|=1. Metric variation precedes FLRW specialization.

## Quantitative contract

In the bulk A2=C²/120-E4/360+R²/72+BoxR/30. The full first variation on CD
is I_R²/72; C vanishes on this background but its second variation need
not. Exactly |ell⁴ I_R²_00|<=6048 and |ell⁴ I_R²_ii|<=4032.

For one independently bounded **constant** scale ratio with
|log(mu/mu0)|<=L, the isolated matter running of cR has magnitude at most
hbar L/[72(4 pi)²]. Its added local source, compared to the positive reference
M²/ell², is bounded by

    168 hbar L / [(4 pi)² (M ell)²]
      <= (7/6) hbar L / (M tau)².

Thus the S5.8 sufficient choice M tau=10^324 gives the conditional bound
(7/6) hbar L 10^-648. No value of L or finite matching coefficient is
inferred from the tree certificate. This is not division by G00, which
vanishes at the bounce. It is not a bound on a solution's displacement.

Finite cR(mu0) is separate matching data; if |cR(mu0)|<=C is independently
established, the corresponding local source ratio is at most
12096 C/(M tau)². No such C is supplied by this audit. A spacetime-dependent
renormalization scale cannot be inserted into these constant-coefficient
variation formulas. Running local terms and the remaining effective action
compensate under a scale change; the isolated increment is not separately
observable.

## Explicit nonclaims

The result does not bound anomaly, massless nonlocal terms, state-dependent
stress, initial/boundary effects, finite matching, other loops, or higher
orders. It neither proves nor disproves physical ghosts, quantum
superluminality, corrected kinetic/gradient positivity or radiative
protection. The fixed spatial momentum band alone does not bound all time
derivatives of a new higher-derivative tensor operator. Corrected coupled
constraints, perturbative order reduction, backreaction and UV matching
require new calculations and error budgets. Full P8 remains open.

The exact certificate verifies algebra and polynomial bounds given the
primary-source heat-kernel theorem. It is not a formal proof of that theorem
or of the uncomputed quantum dynamics. Independent lapse-first and
linearized-Weyl tests supplement the replay.
