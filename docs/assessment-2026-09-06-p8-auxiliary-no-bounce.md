# P8: the full auxiliary parent cannot support the physical bounce

Recorded 2026-09-06, following the
[actual tensor-cone checkpoint](assessment-2026-09-06-p8-auxiliary-tensor-cones.md).
The new result closes the selected exact auxiliary-parent matching route,
including the mixed-sign and asymmetric cases left outside the earlier
cone theorem. It does not close original P8 or change the frozen DHOST rows.

## Main result: use all the background equations

The [S6.15 formulation](../problems/P8/s6/matching/trimetric/cones/global/FORMULATION.md)
keeps two healthy Einstein vierbeine, the actual auxiliary matter vierbein
u, constant real link and endpoint coefficients, and homogeneous matter
acting only on u. At least one link is nonzero. All lapses and scale
factors are smooth, finite and positive on the connected physical domain.
The actual matter null stress is `n_h=epsilon*(rho+p)>=0`; homogeneous
positive-field-metric canonical matter is one sufficient realization.

Let T be u proper time and, for each nonzero link, put

    R_i=a_u/a_i, c_i=R_i*n_i/n_u,
    H_i=(1/n_i)*d(log a_i)/dt, H_u=d(log a_u)/dT.

Full covariant variation, checked again by lapse-retaining action
variation, gives

    3G_i*H_i^2=2p_i*R_i^3+2b_i,
    D_i H_i=(p_i/G_i)*R_i^3*(1/c_i-1).

The Einstein Bianchi identity is not optional additional information. The
interaction density and pressure implied by these equations have balance

    D_i rho_i+3H_i*(rho_i+p_i,eff)
      =6p_i*R_i^3*(R_i*H_u-H_i)/c_i.

It must vanish on an actual solution. Every nonzero link therefore fixes
`H_i=R_i*H_u`, even at zero Hubble rate. Consequently
`R_i'=R_i*(1-c_i)*H_u`. Combining the Einstein acceleration equations
with the actual u null equation yields

    K=sum_nonzero_links G_i/R_i^2 > 0,
    2K*H_u'-H_u*K'=-n_h,
    (H_u/sqrt K)'=-n_h/(2K^(3/2)) <= 0.

No sign assumption on a link, cone inequality, Minkowski vacuum, tail
asymptotic or inverse of the auxiliary tensor Hessian enters this proof.
Disconnected Einstein fields are omitted from K. Smooth positive coframes
make K positive and regular at every finite point; a uniform bound on
both infinite tails is unnecessary.

Because division by positive sqrt K preserves the sign of H_u, a
nonincreasing normalized Hubble rate cannot change from negative to
positive. This also excludes a degenerate crossing with H_u'=0 exactly
at the crossing, or an intervening static interval. At an ordinary bounce
the sharper local statement is simply

    H_u'=-n_h/(2K)<=0 when H_u=0.

An expansion-to-contraction turnaround is not ruled out. Neither is the
previous actual expanding solution; it supplies a nonempty solution
control, not a counterexample to this orientation-sensitive theorem.

## What this changes for matching

The earlier asymmetric low-energy countercontrol remains correct: one
cannot infer a low-energy tensor speed merely from the two full-parent
principal speeds. The new theorem bypasses that issue entirely. There is
no actual parent background with the required physical contraction and
expansion to which a light-only EFT could be faithfully matched.

There is also a finite-interval comparison, conditional on explicitly
identifying the actual physical matter metric, proper clock and endpoints
with those of the old CD target. For

    a_CD(T)=[1+(T/tau)^2]^2,
    H_CD(+-L*tau)=+-4L/[tau*(1+L^2)], L>0,

every regular solution of this exact parent must obey

    max_endpoint |H_u-H_CD| >= 4L/[tau*(1+L^2)].

Otherwise the two parent endpoint signs would be negative and positive,
contradicting monotonicity. At L=1/2 the necessary error is at least
`8/(5tau)`. No field normalization, full scalar action or stationary
vacuum is identified merely by using this physical-geometry comparison.

Likewise the combined equation residual has the conditional lower bound
`2*kappa*a+nu-eta*D` if K>=kappa>0, H_u'>=a>0, n_h>=nu>=0,
|H_u|<=eta and |K'|<=D. A positive result is a necessary cancellation
threshold for a separately specified correction. It is **not** a computed
omitted-operator, loop, state or cutoff error. Adding operators changes
the action and requires its own matching calculation.

## Independent checks of the previously open branches

The complementary cone/tail argument remains informative. For a negative
link, its Friedmann equation already requires b_i>0 and the fixed bound
`R_i^3<=b_i/abs(p_i)`. If that metric cone is inside the matter cone, its
proper-time Hubble rate is nonincreasing. Its scale factor cannot then
grow at both ends as the physical target does. This argument does not
require a vacuum or NEC for that negative link.

A separate genuine mixed-sign flat-vacuum corollary checks the tensor
constraint. Flatness fixes `b_i=-p_i*R_i0^3`; then the same-action lapse
equation implies

    P_i-P_i0=-3G_i*H_i^2/
      [2R_i*R_i0*(R_i^2+R_i*R_i0+R_i0^2)] <= 0,
    P_i=p_i/R_i.

For a mixed-sign vacuum with a finite positive relative Fierz-Pauli mass,
`S0=P_g0+P_f0<0`. Hence `P_g+P_f<=S0<0` throughout every regular
homogeneous solution. The auxiliary tensor denominator cannot cross zero
there. General vacuum ratios require the physical Einstein coefficients
`G_i/R_i0^2`, not the unconverted coefficients. The main background
no-bounce result is stronger and does not need this vacuum corollary.

## Primary-source audit and the remaining research boundary

The common proper-velocity relation has relevant prior art:
Capozziello–Martin-Moruno derive it and simultaneous regular extremality
in equations (15), (18) and (19). Their arbitrary-cosmology reconstruction
allows arbitrary second-sector matter, and their presentation excludes
special parameter branches before (15). It is not a canonical-NEC
bounce witness for the present action.
[arXiv:1211.0214](https://arxiv.org/pdf/1211.0214)

Baccetti–Martin-Moruno–Visser, section 4, establishes the anticorrelation
of paired **interaction** null stresses. Those effective stresses are
not the real matter source. The present full-action, physical-clock
monotonicity and matching bounds are separately derived; that paper is
not cited as already proving this auxiliary theorem or all its branches.
[arXiv:1206.3814v2](https://arxiv.org/pdf/1206.3814)

The mixed-sign vacuum also fails the isolated massive/massless spin-2
eikonal coupling screen: its leading scalar coupling has composite
weights (2,-1). But the source paper permits new physics near the massive
tensor scale instead, assumes the eikonal resummation and a sufficiently
large flat region, and treats the absence of asymptotic time advances as
an additional requirement. This is not an unconditional replacement for
the adopted finite-gravity dispersion contract.
[arXiv:1712.10020, sections 4.5 and 6](https://arxiv.org/pdf/1712.10020)

Falkowski–Isabella's massive-tensor Compton test, equations (4.24)–(4.27),
requires positive scalar coupling in its stated Minkowski EFT. Applying
it to a one-leg decoupling limit still requires an actual amplitude
dictionary and control of the pole/cut and error prescriptions. Its
scalar-elastic equation (4.32) retains massive exchange poles: the sign
of an isolated auxiliary contact is not the complete positivity test.
[arXiv:2001.06800](https://arxiv.org/pdf/2001.06800)

Further vacuum positivity calculations for this unchanged parent cannot
repair its missing bounce background. P8(a)'s photon/global-flat-FLRW
objective and P8(b)'s frozen linear-principal classification remain
complete in their stated scopes. The original UV-matching question stays
open: the new no-go closes an action family, not every possible higher-
operator, noncanonical-source or different-parent realization of the
surviving DHOST rows. No numerical gravitational remainder is invented.

## Verification

The source-frozen S6.15 report has SHA256
`69d032cb10f707410d3d72de3762ef6bdb23b9cc096f11c676fe875917a3f12b`.
All 16 direct source hashes match. It records 58 primary exact residuals,
five independent coefficientwise identities, twelve Fraction/Dual
fixtures and 36 primary/independent comparisons. The separately authored
27-test audit is included in the frozen inventory.

All 65 ordinary tests passed in 117.91 seconds; the standalone read-only
certificate replay and Ruff passed. The written monotonicity argument
supplies the continuous-time implication; these arithmetic checks are
supporting evidence, not a proof-assistant formalization.
