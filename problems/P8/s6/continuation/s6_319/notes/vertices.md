# General vertex and hard-path estimates

## All-valence vertices and hard-path propagators

Every light scalar momentum on the core has absolute components<3,
hence Euclidean norm<6. Every heavy path momentum has components<5
and Euclidean norm<10. These follow from external COM bounds plus
a subset of total radiation, independently of N.

S315 gives the r-metric scalar vertex bound
 r!2^r[binom(r+2,2)||p||||q||+(r+1)m2]*product||A||.
For r>=1 use binom(r+2,2)<=(r+1)^2 and (r+1)^2<=4^r.
For light scalars the bracket is<=37(r+1)^2, and37*8^r<=1024^r.
For heavy scalars, n>=128 and100/n<=1 give a bracket bounded by
2n(r+1)^2, and2*8^r<=1024^r. Thus the vertex envelopes are
 r!1024^r, or n*r!1024^r, with canonical powers restored.
Density vertices obey r!2^r(r+1)<=r!4^r for r>=0.

Every extra H metric vertex is paired with one heavy inverse<2/n.
The n in its vertex cancels, leaving one overall heavy inverse2/n
for the path, irrespective of its length. No high-mass pole is crossed.

A hard h path has one fixed underlying Born channel. On a mixed path,
S312/S315 give
 |D_next|>(tau+W^2)/300000,
 max_component(momentum)<=sqrt(tau)+4W.
The maximum includes the two hard legs and every attached block
momentum. Hence each squared Euclidean momentum norm is at most
4*(sqrt(tau)+4W)^2<=128*(tau+W^2).
Their ratio to the NEXT hard propagator is<38400000.
On a timelike path max_component<=4 and D_next>=45/8, so the
same common bound is conservative.

The four-dimensional trace-reversal map on symmetric covariant tensors
is an exact Euclidean Frobenius isometry: it reflects the normalized
eta direction. Raising/lowering indices with eta is also an isometry.
Consequently an EH vertex with two hard legs and r>=1 blocks, followed
by the next hard inverse, has normalized multilinear operator bound
 38400000*(r+2)!*32^(r+2).
This uses the full vertex r+2 bound and no conservation of a collapsed
hard source. It holds for complex tensors using Hermitian norms.

The one initial hard inverse is bounded by300000/(delta+W^2).
For timelike channels the same estimate is far looser than16/45.
After division by AG>8/(K*delta), its relative contribution is at most
300000/8, including the single extra hard K^(-1). Every additional
EH step contributes only the canonical powers of its soft attachments.
There is no per-step accumulation of forward poles.
