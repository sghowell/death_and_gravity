# P8: full finite-kappa analytic-target classical match

Original P8 remains OPEN. This follows the
[complete two-loop low-energy cut](assessment-2026-09-10-p8-complete-two-loop-elastic-cut.md).
No user intervention is required for the current research.

## From a quartic target to the complete analytic target

[S6.162](../problems/P8/s6/continuation/s6_162/FORMULATION.md)
strengthens the earlier restricted classical stationary-action
match. It now compares against the FULL finite-kappa S6.109
rank-regular analytic flat scalar target, not only its
mass-one quadratic and quartic germ.

The target retains its actual rational-step switch,
n=1024, kappa=10^800 and lambda=10^-600. The canonical
conversion is u=Phi/sqrt(kappa), X=Y/kappa in the fixed
light-mass units. In particular a dimensionless unit-clock
gradient is not a unit canonical gradient.

The flat canonical density is

    L=kappa F+(A3 L3+A4 L4)/kappa+A5 L5/kappa^2.

The curvature term vanishes only on the stipulated fixed
flat metric. Its metric variation is not set to zero.

The exact sixth/eighth field-degree germs are extracted
from the literal target. The sixth density is

    -(2681 Phi^2+853Y)(Phi^4+nY^2)/(200 kappa^2)
    +3n(2Phi^2+Y)(L3-L4)/kappa^2.

It is nonzero and is not absorbed into a quartic identity.
The report and written proof retain the complete eighth
density, including all lower-scalar and DHOST groups.

## All remaining field degrees and a finite action bound

For the 21 ordered canonical scalar/first/second jets,
the unit-cube invariant caps are

    |Y|<=4, |L3|,|L4|<=64, |L5|<=256.

Exact coefficient majorants give

    E6=2.858465025e-1594,
    0<E8<3e-2195.

For a common complex field-amplitude parameter t, choose
|t|=rho=10^300. Then |u|<=1e-100 and |X|<=4e-200.
The full rational switch is bounded, not discarded.
Every target denominator has a strict gap, including
|R-1|<1/2, and the complete canonical density is
holomorphic on a neighborhood of the closed circle
with modulus below M=10^811.

Parity makes the density even in t. Cauchy's estimate
therefore bounds every unretained degree at least ten by

    E10=M/[rho^10(1-rho^-2)] <2e-2189.

This is a field-amplitude estimate, not a momentum
cutoff or quantum-loop remainder.

At a physical point let W be the maximum modulus of
the 21 canonical jets. Rescaling the jets by W gives

    |L_full-L2-L4|<=E6 W^6+E8 W^8+E10 W^10
                   <=(E6+E8+E10)W^2 for W<=1.

Use the SAME real Schwartz/Fourier unit class of S6.111,
with Fourier L1 norm at most one and U=||Psi||_2.
Parseval gives integral W^2<=21 U^2, so the full
higher-field target action allowance is

    21(E6+E8+E10)U^2 <1e-1591 U^2.

The coefficient is approximately 6.0027765525e-1593.
No constant pointwise bound is integrated over infinite
spacetime.

Adding this allowance to the previous field-map and
unexpanded heavy-resolvent error proves

    |S_polynomial[F(Psi),H_stationary]-S_full_target[Psi]|
      <1e-800 ||Psi||_2^2.

The mapped field/source support radii remain 3 and 6;
the centered operator radius remains 38. The final
strict inequality uses exact rationals.

## Independent verification

The native report pins 18 source/proof/test files and
20 fields: 40 named identities, 40 scalar entries,
31 proof gates, 9 controls and 174 rejected inputs.
Final private science passed 294 tests in 6.40 seconds;
fresh repository science passed 294 in 6.61 seconds.
Ordinary replay passed 319 tests in 2342.01 seconds.
Independent native CLI replay passed. The complete
captured P8 snapshot passed 26229 tests in 3345.01
seconds with final exit code 0.

All 577 captured test files were present and unchanged.
Path-list SHA-256:

    93dfbe1dd182d4480c0b0173efb99d83758261c658ed686856e3f4efd6631212

The full-run adapter passed 128 original tuple comparisons.
Counters were 34099 domain fallbacks, 7340 exact descents
and 88 mixed fallbacks. This snapshot predates S6.163.
Only full regression uses the separately audited exact
GCD adapter. Native, direct science, ordinary and CLI
retain unmodified SymPy with interpreter-only allowances.

The 75661-character native report was transferred
losslessly in seven chunks, and all eighteen source
hashes were independently verified. Report SHA-256:

    50d2439dea8d9b3889f4875679c05f27d2d4e9e81ef88d33bc404cb4671e3cfd

Independent tests extract coefficients of the full
rational/exponential functions on complex circles,
evaluate the actual calibrated target at 2450-digit
precision, check literal Lorentz contractions and
test the Cauchy/integrated-error arguments.

Cleanup preceded final private verification and freeze.
No frozen scientific source or report changed. Written
analytic arguments and exact replay are not formalization
or independent peer review. Exact 22-file staging excludes
later continuations and unrelated P4/P9 changes.

## Remaining work

The separately named SAT8 protected-mass profile and a
polynomial clock-transparent source coordinate have passed
native and fresh direct science; their independent
regressions are running. Neither is counted as published
in this checkpoint. A free flat-clock Dirac production
estimate is being developed with an explicit in/out state
and transition-series remainder.

The classical full-target estimate does not imply quantum
target matching, metric/stress control, a global inverse
or a common rolling parent. Physical finite-EFT truncation,
global V contours/cuts, finite-gravity G and common-parent
B remain open. No all-orders UV requirement is added.
Scoped P8(a) and A.20-A.23 are unchanged; original P8(b)
and original P8 remain OPEN.
