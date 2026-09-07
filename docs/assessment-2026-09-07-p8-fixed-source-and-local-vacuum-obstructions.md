# P8: fixed-source and local-vacuum limits of the current parent

This continues the [causal-response and reciprocal-boundary checkpoint](assessment-2026-09-07-p8-causal-response-and-global-boundary.md).
The next two results test two tempting inferences about the same local
parent: that its large center algebraic mass guarantees a unique
fixed-source limiting response, and that its regular bounce automatically
supplies a vacuum for the adopted UV test. Both inferences fail on the
explicit domains considered here. Neither failure classifies every EFT
or closes original P8.

## A fixed smooth source can retain a phase-dependent response

The [fixed-pulse theorem](../problems/P8/s6/matching/variable/response/forced/phase/FORMULATION.md)
starts with all four physical tensor data zero at u=-r. Fix
`0<r<=1/1000`, spatial momentum `0<=K<=4`, and a real source sigma that
does not depend on delta. It lies between zero and one, is supported in
`[-r,-r/2]`, and differs from the unit rectangular pulse by an L1 deficit
at most r/100. This class includes explicitly constructed smooth pulses
with compact support strictly inside that interval.

For the unchanged full coupled action and actual physical g observable,
the endpoint response remains bounded as delta tends to zero, but

```
limsup gamma_g,delta(r)-liminf gamma_g,delta(r)>3*r^2/20.
```

The source, initial-data prescription, spatial momentum, measurement
time and physical metric are all held fixed. This is not a numerical
scan of oscillations or a comparison of different source preparations.
Scaling the source by any fixed positive amplitude scales the separation
by the same factor, within linear response.

The exact finite-clock heavy kernel gives an oscillatory term with phase
`(sqrt(39)/2)*log(32*r^2/delta)`. Its amplitude, after retaining the actual
source and observable weights, is greater than `2*r^2/25`. The full
coupled feedback and bounded pole remainder change the endpoint by less
than `1200*r^4`. These constants leave a strict separation between two
sequences of delta values. The actual residuals need not converge on
either sequence for the limsup/liminf conclusion to follow.

The Gauss connection is used between its zero and one arguments, with
the correct Euler transformation for the second left-moving factor.
An all-degree coefficient-ratio bound controls the tails uniformly on
the fixed source interval. Physical-time normalization, both clock
phases and the source moment are retained. The result survives the
specified smooth turn-on and turn-off; it does not rely on a distributional
impulse or a discontinuous source.

This excludes convergence to a delta-independent response for this pulse
class. It does not exclude an explicitly delta-dependent, phase-retaining
description. Nor does it establish that every low-energy EFT fails:
these finite compact pulses do not have a parametrically small version
of the previously defined RMS/stiffness-proxy hierarchy. The RMS claim
applies to H1_0 members of the class; it need not be defined for general
bounded pulses. Long or approximately bandlimited preparations and
selected correlated initial data remain distinct questions.

## The same local action has no proportional constant-clock vacuum

The [local vacuum theorem](../problems/P8/s6/matching/variable/vacuum/FORMULATION.md)
tests the full stationary equations, not an instantaneous matrix on the
bounce. Keep the original positive Einstein terms, both canonical
scalars and all original beta functions. Let both scalars be constant,
`g=eta` and `f=r_vac^2*eta`, with every positive proportional ratio
r_vac considered independently of the action parameter c.

The clock value is allowed anywhere in the original local field interval,
defined by `phi=M*integral_0^u sqrt(kbar)du`, `|u|<=1/10`. Here u labels
a coefficient value; it is not a time coordinate in the proposed vacuum.
The quantities used in the rolling reconstruction are not incorrectly
imposed as nonzero kinetic stresses in the stationary configuration.

Literal metric and clock variations give, for the actual beta2=beta3=0
profiles, the three simultaneous conditions

```
b0+3*r_vac*b1=0,
b1+r_vac^3*b4=0,
[2M/(tau^2*sqrt(kbar))]*(b0,u+4*r_vac*b1,u+r_vac^4*b4,u)=0.
```

Each metric's spatial equations are also checked. The clock equation's
four-dimensional binomial weights differ from the three-dimensional
weights in the lapse equations; omitting that distinction would change
the problem.

The second metric equation has a unique positive root throughout the
stated field interval. At that root, the physical-g residual satisfies
the continuous bounds

```
U_vac=b0+3*r_vac*b1 < -123/25       for c=1,
U_vac=b0+3*r_vac*b1 < -459/800      for 2<c<=4.
```

The physical density is `2M^2*U_vac/tau^2`, whereas the flat Einstein
tensor and constant-scalar kinetic densities vanish. Thus the metric
equations cannot be solved simultaneously, regardless of the clock
equation. At the center, r_vac=2 satisfies the second-metric and clock
equations but leaves `U_vac=2-16/c<0`; that is not a vacuum.

The proof uses exact rational interval bounds, a residual identity and
a continuous lapse polynomial. It is not an exclusion inferred from a
root-finding grid. The density gaps are conditional on imposing the
second-metric equation exactly, not a generic error budget for arbitrary
approximately stationary configurations.

This fails the adopted V test's algebraic vacuum prerequisite only within
the specified local field domain. It does not exclude a separately named
smooth extension with a vacuum at another field value. Adding a constant
to b0 on the original interval can repair a center vacuum equation, but
then changes the original bounce equations; the certificate includes a
control exposing that action change. No vacuum spectrum or positivity
amplitude is inferred from the absence theorem.

## Remaining research and verification

The current parent now has several precisely located limitations, rather
than a claimed universal decoupling or UV verdict. A new vacuum extension
would still need its complete canonical spectrum and the same-parent B
matching controls; a new source class would still need a quantitative
state and omission estimate. Neither a mathematical extension nor a
selected light subspace alone supplies the original C/D operator
dictionary, cutoff, omitted-loop control or finite-gravity dispersion
remainder. The adopted contract does not require a global physical-time
history connecting a vacuum to the bounce.

The completed scoped photon objective and original 32 linear row verdicts
are unchanged. Original P8 remains open, and no new user decision or
permission is required to continue the research.

Both new proofs have passed independently authored scientific audits.
The phase audit passes 22 tests and includes literal source normalization,
Euler/Gauss tails, clock factors, continuous moment bounds and limsup
logic. The vacuum audit passes 15 tests and independently varies all
four g/f coframes before imposing proportionality, checks the full clock
equation and reconstructs the continuous obstruction bounds.

Both source manifests are frozen and all recorded source hashes have
been checked against the files. The phase report hashes to
`8e5185adf4023d7a656a9b99b084b7a3adf76a176c9a0eaaaa692f187c6b6660`
and covers 12 sources, eight exact identities, 11 continuous margins,
108 independent rational comparisons and 20 invalid-input controls.
Its 64 ordinary tests passed in 162.40 seconds, and a separate fresh
ordinary `--check` replay passed. The vacuum report hashes to
`ad04e78149b4c21fe01c4b093fa4b718d0b85bb8c5634d1d6092b0b380b342a1`
and covers 14 sources, 23 identities, 171 independent comparisons,
three Bernstein coefficients and 17 invalid-input controls. Its 58
ordinary tests passed in 147.12 seconds; its separate fresh ordinary
`--check` replay passed in 128.41 seconds.

The expanded P8 regression passed all 2,834 tests in 613.52 seconds,
with no exclusions. The exact GCD regression adapter passed its 128
self-checks and recorded 6,509 exact descents and 4,527 original-domain
fallbacks. Ordinary certificate replays do not use that adapter. Both
Python's hash seed and SymPy's separate polynomial-factorization random
seed were fixed to zero. The host's problematic faulthandler timer was
disabled, without skipping mathematical assertions. Frozen ancestors
and unrelated P4/P9 changes remain untouched.
