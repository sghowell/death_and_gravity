# P8(a) A.14: a longer solution with the identical preparation and state

This checkpoint proves an actual extension of A.11, not merely a small
frozen-operator bound. It uses A.13's full causal resolvent, including its
growing pole; it imposes no future boundary or no-runaway condition.

## Unchanged physical inputs and exact source

Keep the one free massless minimally coupled scalar, original positive
quasifree Hadamard vacuum transported from the unchanged radiation past,
fixed ordinary radiation, zero cosmological term, no additional explicit
curvature-squared coupling and the named prescription

    lambda=2sqrt(2) A eta_star^2, gamma=0,
    epsilon=1, delta=10^-14,
    kappa hbar=2880pi^2 delta A^2 eta_star^4.

Use the exact conformal coordinate x=eta/eta_star-1, a_phys=A eta_star a,
u=-a''/a and h=a'/a. Let x0 be the old plateau point y0=5/2 and t=x-x0.
The baseline (barred) geometry and actual state are the same old prepared
ones, continued by their original formulas on the plateau; no inverse
clock expansion or replacement quantum covariance is used.

The old preparation duration remains **L0=10^-10**. Use exactly the same
chosen smooth cutoff chi_L0 as A.11, equal to one on [0,L0/4] and zero
from L0/2 onward. Set c=chi_L0*cbar, with cbar the same fixed weighted
actual old density defect. Do not rescale this cutoff to the new interval.
The exactly conserved source is still

    rho_ext=F c/a^4,
    p_ext=F[c/(3a^4)-c'/(3h a^4)], F=hbar/(pi^2 A^4 eta_star^8).

It is unchanged in the original past and zero after the original cutoff.
As in A.11, it is not globally compact in the off-shell radiation past.
The ordinary radiation normalization and quantum integration constant are
not adjusted. Initial metric, Wick and energy-constraint data are identical.

## The longer actual solution

On the new future interval [0,L], with **L=10^-6**, put

    X=u'-ubar', W=I X, sigma=2*10^6,
    ||X||_sigma=sup_[0,L] exp(-sigma t)|X(t)|,
    X(0)=0, ||X||_sigma<=r=2*10^-7.

The full resummed map is a strict contraction and preserves this ball.
Its fixed point is smooth and defines the actual transported Hadamard
state. The same conserved full SEE holds during preparation, and the
**full unforced SEE holds on x0+L0/2<x<x0+L**. This is a 10^4 increase
in the constructed total future length, with no increase in source duration.

Every ball input satisfies pointwise barriers

    2<=a<=3, 1/3<=h<=1/2, |u|<3*10^-12,
    |u'|<=10^-5 on the whole original history of duration<=3,
    |P=a^2q'|<10^-5, |q|<2*10^-8.

The exact weighted full-map Lipschitz constant is below 9*10^-6 before
the inverse; the resulting contraction is below 3*10^-6. The center
image is below 9*10^-8, and the fixed point has

    ||X*||_sigma<9*10^-8, ||X*||_infinity<81*10^-8.

The exponential factor exp(sigma L)=exp(2)<9 is included in these
pointwise conclusions. A weighted-ball bound alone is not substituted
for them. The new solution lies in the old A.11 uniqueness ball on the
old interval, solves the identical original map there, and hence agrees
with A.11 exactly. It is not a different preparation branch.

## Verification and limits

The certificate pins/replays A.13 and A.11, the full weighted constants,
pointwise barriers, complete source terms, independent Fraction arithmetic,
and the written smooth-state/constraint proof. The infinite-frequency
response, functional analysis and Hadamard implications remain written
proofs, not conclusions from finite mode sampling.

No longer-time stability, growing-pole suppression, all-order jet-norm cap,
physical validity of the switch, realistic-field theorem or cosmological
incompleteness is established. The A.12 QSEI coefficient/reference credit
is not automatically transferred to this larger domain. A fresh QSEI and
focusing-duration calibration is the next task. P8(a) and P8 remain open.
