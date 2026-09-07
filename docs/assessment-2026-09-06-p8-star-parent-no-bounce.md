# P8: closing the constant star-parent bounce route

Date: 2026-09-06. This extends the
[auxiliary-parent no-bounce assessment](assessment-2026-09-06-p8-auxiliary-no-bounce.md)
to a larger, explicitly specified parent-action class. The completed
linear classification and the scoped photon result are not reopened.
Original P8(b)'s common-parent and UV-applicability problem remains open.

## Result and action boundary

Changing constant interaction coefficients, including their signs, or
adding a healthy central Einstein term does not repair the bounce within
this family. The new [S6.16.STAR formulation](../problems/P8/s6/matching/star/FORMULATION.md)
uses a finite star of pairwise square-root interactions,

```
S = -G_u/2 integral sqrt|h| R_B[h]
    -sum_i G_i/2 integral sqrt|g_i| R_B[g_i]
    -2 sum_i integral sqrt|g_i| sum_(n=0)^4 beta_in e_n(sqrt(g_i^-1 h))
    + S_m[h].
```

All coefficients are constant and real, `G_i>0`, and `G_u>=0`.
The only actual matter is homogeneous/isotropic central matter with
`n_h=rho_m+p_m>=0`. A separate central cosmological term is allowed.
The metrics share a smooth, positive, regular spatially flat FLRW/root
chart. There is either a positive central kinetic coefficient or at least
one genuine interaction link. No tensor-cone, flat-vacuum, heavy-gap or
auxiliary-Hessian condition is assumed.

On every connected regular interval the theorem is

```
H_u(T0)<=0  implies  H_u(T)<=0 for every later T.
```

Thus there is no contraction-to-expansion transition, including a
degenerate transition or a zero plateau. This is not a theorem of
geodesic incompleteness, nor an assertion of perturbative stability.

## Why the branch analysis matters

Write `R_i=a_u/a_i`, `N_i=n_u/n_i`, `c_i=R_i/N_i`, with physical time
`dT=n_u dt`. Literal independent lapse and scale variations give

```
J_i(R)=beta_i1+2 beta_i2 R+beta_i3 R^2,
6 N_i J_i(R_i)(R_i H_u-H_i)=0.
```

The Bianchi equation must not be divided by `J_i` everywhere.

| Leaf type | Actual background consequence | Kinetic combination |
| --- | --- | --- |
| Dynamic, `J_i!=0` | `H_i=R_i H_u` and `R_i'=R_i(1-c_i)H_u` | Includes `G_i/R_i^2` |
| Genuine persistent root, `J_i=0` | Constant `R_i,H_i`; `H_u'=(c_i'/c_i)H_u` | Does not include this leaf |
| Endpoint-only, `beta_i1..3=0` | Separate cosmological terms; no clock-lock constraint | Omitted |

For a fixed dynamic set D with positive
`K_D=G_u+sum_D G_i/R_i^2`, the full Einstein equations imply
`(H_u/sqrt(K_D))'=-n_h/(2 K_D^(3/2))`.
There is not generally one exact all-leaf K across algebraic interiors.

The [written proof](../problems/P8/s6/matching/star/notes/proof.md)
instead handles arbitrary branch zero sets. The clock-lock function and
its derivative vanish on dynamic open sets and, by continuity, at their
boundaries. On algebraic interiors the direct sign relation above applies.
On any compact regular interval, positive smooth ratios and lapses give
a finite logarithmic-rate bound C. Both cases imply `H_u'<=C H_u`
where `H_u>=0`. The absolutely continuous positive part
`Y=max(H_u,0)` obeys `Y'<=C Y`; an integrating factor proves no crossing.
No finite-switch assumption or numerical sampling of possible histories
is used. No continuation through a singular lapse or root is asserted.

## Actual controls, not just formal branch labels

An exact algebraic solution has

```
G_u=G_i=1, beta=(0,1,-1/2,0,1/2), T>0,
a_u=a_i=T^(1/3), n_u=1, n_i=1/(3T),
phi=sqrt(2/3) log T, V=0.
```

Both metrics' independent Einstein equations and the scalar equation
hold. The leaf has `H_i=1`, the physical metric `H_u=1/(3T)`;
the interaction null densities vanish. Using the incorrect all-leaf
`K=2` produces the nonzero residual `-2/(3T^2)`. The correct dynamic
coefficient is only `K_D=G_u=1`. This expanding background neither
contradicts the sign theorem nor proves this algebraic branch healthy.

The nontriviality hypothesis is necessary. With `G_u=0`, all interactions
zero, a Minkowski leaf and zero matter, the center has no equation of
motion. Its scale factor can be `1+T^2`. Independent literal potential
variations and Einstein/source residuals verify this undetermined
exception. It is not a gravitational bounce model; adding `G_u>0`
restores a central Einstein constraint that excludes it.

The original auxiliary action is a literal subfamily:
`beta_i3=p_i`, `beta_i0=b_i`, other leaf betas zero, `G_u=0`, and the
central potential `-2B sqrt|h|`. Then `J_i=p_i R_i^2`, so every nonzero
original link is genuine and has no positive algebraic root.

## What is closed and what remains

This closes the regular constant-star, single-central-NEC parent route
to an actual flat physical bounce. It is stronger than excluding a
particular vacuum, cone sign or inverse-Hessian patch. Matching a formal
coefficient at a point does not evade the obstruction when matching also
claims an actual common-parent background and the stated physical-clock
dictionary. The auxiliary child's explicit two-endpoint CD error bound
provides one quantitative instance of this distinction.

This does not exclude arbitrary DHOST rows or all local UV completions.
Variable couplings, derivative interactions, extra matter sectors,
different interaction graphs or an actual non-NEC quantum source require
their own equations and matching/error analysis. Some such changes would
also change the adopted problem contract; none is silently accepted as a
solution here. A tree-level zero, a generic higher-operator threshold or
an uncomputed loop correction is not a UV verdict.

The [primary-source audit](../problems/P8/s6/matching/star/notes/sources.md)
separates the literature's interaction-NEC reciprocity and regular-clock
results from this action's actual-matter assumption and its all-branch
proof. Exact arithmetic and independent coframe variations support the
written argument; it is not proof-assistant formalized.

## Verification

The source-frozen S6.16 report has SHA256
`54f23dfd05d7580969550d06ac4a72adc5faa61988be9e7c95eeb199fb546353`.
All 16 source hashes match. It contains 26 primary exact residuals,
ten independent coefficientwise identities, 24 literal coframe-jet
directions, fourteen Fraction fixtures and 34 primary/independent
comparisons. The separately authored 28-test audit is source-pinned.

All 77 ordinary tests passed in 119.83 seconds; the standalone read-only
certificate replay and Ruff passed. A combined fresh run of the two
new independently authored audit files passed all 55 tests in 1.21 seconds.
No ancestor certificate was regenerated and no unrelated P4/P9 work
was included in these gates.

The combined P8 regression passed **2,040 tests in 424.31 seconds** using
the explicitly opt-in [exact regression runner](p8-exact-regression-runner.md):
128 comparisons with original Gaussian-domain GCD tuples passed at startup,
followed by 6,509 exact descents and 2,503 domain fallbacks. This is the
adapted exact full-suite run, not an ordinary unmodified SymPy run.
All ten ordinary runner tests also passed. One unchanged long-running
symbolic test emitted the configured 60-second traceback diagnostic and
then completed successfully; it was not counted as a stopped-run pass.
