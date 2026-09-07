# P8: uniform scalar principal control, with a wider full-parent cone

This continues the [prepared-tensor checkpoint](assessment-2026-09-07-p8-prepared-tensor-sector.md).
The unchanged local variable-beta parent now has a uniform scalar
high-frequency theorem through the bounce, not just regular constraints
or a frozen center matrix. It does not yet have a controlled general
light-only EFT, a global scalar-health theorem, an original C/D operator
match, or the adopted vacuum/finite-gravity UV verdict. P8 remains open.

## The scalar issue that is now resolved

S6.22 proved a regular six-state scalar Hamiltonian but showed that one
physical observable chart changed momentum degree at the center. Its
naive frozen-center value therefore did not identify a uniform cone.
The [S6.24 formulation](../problems/P8/s6/matching/variable/perturbations/microlocal/FORMULATION.md)
and [complete proof](../problems/P8/s6/matching/variable/perturbations/microlocal/notes/review.md)
resolve that issue for the specified local constant-lapse family.

The new physical variables are the relative spatial scalar `pi=-e/K`,
the g-clock displacement xi and the actual Bardeen matter perturbation
chi_B. Both the literal unreduced action and the exact physical Cauchy
map retain all time derivatives before evaluation at the bounce. A
boundary coefficient that vanishes at the center nevertheless has a
nonzero derivative there. At c=4 that derivative changes the helicity
gradient from the incorrect 448 to 800. The correct kinetic coefficient
is 480, giving squared speed `5/3`, not the discarded `14/15` or the old
nonuniform chart's value 5.

For each fixed `0<delta_min<=2`, on `|u=T/tau|<=1/10` and
`2+delta_min<=c<=4`, the exact Cauchy map is uniformly regular for

```
K=(tau*k_com)^2 >= 2000000/delta_min^2.
```

This deliberately loose bound follows from the exact cubic determinant
and continuous rational coefficient bounds. It is a sufficient
mathematical symbol threshold, not a physical cutoff or a demonstrated
range of EFT validity. Exact polynomial identities show that the full
physical equation is `q''=-K K_s^-1 G_s q+O(1)q+O(1)q'`, uniformly on
each fixed-delta box. The proof controls the compactified `1/K=0` limit,
rather than extrapolating from finite samples or a generic degree count.

After the time-dependent weight `Pi=(a^2*y*P/2)pi`, both principal matrices
obey `(1/4)I < Kbar,Gbar <64I`. Four continuous Bernstein certificates
contain 611 strictly positive coefficients. The normalization's first
and second time derivatives are retained in the exact evolution.
The resulting Fourier energy obeys

```
E(u) <= exp(C_delta*abs(u-u0))*E(u0),
```

uniformly in all K above the displayed threshold. C_delta is a finite
explicit-formula compact supremum, not a reported numerical optimum.
The conserved symplectic form independently identifies the positive
principal kinetic weight. This is stronger than choosing an artificial
positive norm for an otherwise unverified system.

## What this does and does not say about causality

The three center squared scalar speeds in the actual physical g frame
are `1,1,(9c^2-22c+144)/120`. The last is strictly greater than `17/15`
for every `2<c<=4`, and equals `5/3` at c=4. Thus this full-parent scalar
principal sector is positive but not fully subluminal.

The theorem holds at fixed positive c-2. Although the weighted principal
matrices extend smoothly to the c=2 coefficient edge, the actual parent
and its time connections do not. The sufficient threshold grows as
`delta_min^-2`, whereas the center algebraic tensor mass squared grows
only as `(c-2)^-1`. This high-frequency proof cannot be imported as a
subthreshold heavy-mode elimination theorem.

The energy uses explicit, fixed spatial derivative weights. It is not a
uniform bound on every raw metric amplitude, an instantaneous ground
state, all-band/nonlinear stability or a prescribed-source response.
It also does not transfer to the different, variable-lapse G1 action.

## Next research gates

The selected analytic tensor sector remains distinct from the general
retarded response to physical matter. The next tensor calculation is a
causal, source-normalized elimination with explicit state, time-window
and omission bounds. A nonlocal response formula must not be relabeled
as a local low-energy EFT without a justified derivative hierarchy.

Global-parent work must separately retain its full metric domain and
actual scalar/vector dynamics. Neither physical-g completeness alone
nor an additional metric's incomplete chart decides the original EFT
classification. The common-parent C/D dictionary, controlled spectrum
and cutoff, omitted operators/loops and adopted V/G UV hypotheses remain
required where applicable. The completed scoped photon objective and
32 original linear row verdicts are unchanged. No new user choice or
permission currently blocks the remaining research.

## Verification record

The frozen [report](../problems/P8/s6/matching/variable/perturbations/microlocal/certificates/variable-scalar-microlocal.json)
has SHA-256
`56392e03143c4658170826565519c7ed3ff9a90ebc2140cfa8c1ec8cdb52e65a`.
It pins 21 sources and S6.22's recursively verified ancestry, with 188
exact action/completion identities, 36 exact energy identities, 611
continuous positive Bernstein coefficients, 1,373 independent Fraction
comparisons and 13 rejected-input controls. The separately authored root
audit passes 26 tests, including all continuous matrix bounds, the
general-center determinant and full normalized second-jet Cauchy map.

All 59 final ordinary tests passed in 224.83 seconds. A separate fresh
ordinary CLI rebuilt the report and exited successfully with Python's
hash seed and SymPy's independent random seed both fixed to zero. All
21 source hashes remained unchanged. The expanded regression passed
2,592 tests in 586.23 seconds using the documented exact-arithmetic
adapter: 128 startup self-checks, 6,509 exact real descents and 4,519
out-of-domain fallbacks. Ruff and whitespace checks passed. Interrupted
host-diagnostic runs are not counted as successful verifications.
Unfinished reciprocal-geometry and forced-response children were excluded
from this checkpoint manifest and regression; unrelated P4/P9 work was
preserved.
