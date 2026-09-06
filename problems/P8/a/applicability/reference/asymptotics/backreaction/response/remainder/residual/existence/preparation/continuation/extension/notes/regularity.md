# Same-data smoothness, actual-state identification and old uniqueness

The C0 weighted gate alone is not the final physical theorem. The
remaining steps are justified here using the already proved A.11
flat-start response lemmas and A.13's complete causal kernel Kf.

## 1. No initial endpoint or state change

Start Picard iteration at X0=0. Every iterate is zero on the **same
original** neighborhood [0,L0/4], since c=cbar there and the full metric,
auxiliary and mode histories then coincide. Causality preserves this
property. All iterates are smooth and agree with the original past with
every jet unchanged. This uses the original cutoff, not a new smooth
join at a moved endpoint and not a reset of covariance on a new metric.

Kf is locally L1 and has weighted norm <=k=240/769. For smooth input f
zero on that fixed neighborhood, differentiation can be placed on f:

    (Kf*f)^(m)=Kf*f^(m),
    ||Kf*f^(m)||_sigma<=k ||f^(m)||_sigma.

No derivative of the singular kernel, initial boundary term or generic
C1-endomorphism assertion is used. Smooth extensions past L do not
affect this retarded statement on the interval.

## 2. Uniform highest-jet block in the weighted norm

The actual A.11/A.10 polarized Dyson argument gives, for each finite m,

    partial_t^m R'[u]=DR[u](u^(m+1))+lower jets,
    ||DR[u]Z||<= [z+18z^2 exp(2z)]||Z||, z=MT^3.

The quadratic term is its convergent polarized logarithmic triangle;
it is not replaced by a divergent absolute ultraviolet integral. Higher
orders use the summable factorial Dyson majorants. The operator is causal.
For directions Z zero in the common past, apply its bound on each prefix
and then multiply by exp(-sigma t). This gives the identical coefficient
in the weighted norm. Operator-norm continuity in u follows from A.11's
explicit second-polarization bound; equivalence of the finite-interval
weighted and unweighted norms supplies any needed finite factor E.

The differentiated new Picard map has

    X_(n+1)^(m)=A_n X_n^(m)+B_(m,n),
    A_n Z=-Kf*[(d(a_n)-d(af))Z+DR[u_n]Z].

All other terms in the exact forced remainder involve lower X jets:
u gains one derivative, h gains two, a gains three, P gains two and q
gains three. Source derivatives are fixed smooth functions, although
their high-order bounds may be large. The Einstein term is either in
the fixed convolution inverse or in an integrated rolling remainder;
it does not multiply the highest X^(m) in B_(m,n).

Consequently, independently of m,

    ||A_n||_sigma<=k[Df+z+36z^2]<10^-4<1.

The lower remainders are bounded when lower jets are bounded and converge
when lower jets converge, by the same finite ODE/product and polarized
Dyson arguments as A.11. For m=1,2 this is read directly from the
integral equations; no negative derivative orders are asserted.

## 3. Smooth convergence on the longer fixed slab

C0 convergence is already proved. Inductively assume convergence through
order m-1. The last display bounds the mth derivatives uniformly by a
strict scalar recurrence with bounded lower remainder. Differences of
two late iterates obey that same contraction plus an error tending to
zero, using operator-norm continuity and lower-jet convergence. Taking
the limiting supremum over late pairs proves that the mth derivatives
are Cauchy. The norms are equivalent on this finite interval, so uniform
derivative convergence identifies their limit. This proves C-infinity
convergence for every finite order on **the same [0,L]**, not shrinking
intervals or a quantitative bound uniform over all derivative orders.

## 4. Full stress, constraint and identity with the old solution

The smooth positive metric transports the original Hadamard vacuum
by the same Klein--Gordon evolution. A.11's actual-mode/Hadamard
identification applies without changing its hypotheses: the past is
identical and all jets match. The now smooth actual Wick functional
coincides with the auxiliary one because their derivatives agree at
the fixed point and their initial values coincide. The exact anomaly
and named finite scheme therefore yield the actual stress tensor.

The source remains exactly conserved. Its same initial value supplies
the same full energy constraint; the inherited exact conservation
identity propagates that constraint. Thus both density and trace SEE
hold pointwise, and after t=L0/2 the external source is identically zero.
This is an actual longer source-free solution with the fixed ordinary
radiation normalization, not a frozen-background residual statement.

The new solution solves the identical original A.11 fixed-point equation
on [0,L0], since resummation is an exact causal identity and the source,
state and data have not changed. Its global pointwise bound is already
<9*10^-7<10^-6, so it lies in the old uniqueness ball on that interval.
A.11's uniqueness proves exact agreement there. No extra homogeneous
pole amplitude or no-growing future condition is selected.

This does not transfer A.12's reference positivity or its QSEI constant
onto the longer interval. Nor does it give analytic jet bounds, a large
time stable manifold, field fluctuations or cosmological incompleteness.
