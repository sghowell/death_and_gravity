# Spatial Sobolev and full source derivative bounds

Let I have length at most one. The real scalar history Phi is
smooth, spatially Schwartz and zero on an initial neighborhood.
All its coordinate jets through order three, including mixed
time/spatial jets, have pointwise absolute value<=1.
For any smooth function f write

    U_f(t)=max_(|alpha|<=3) ||partial^alpha f(t,.)||L2(R3),
    ||f||J3=||U_f||L2(I).

No Fourier band limit is assumed, and no compact time support
inside the final observation endpoint is required for Phi.
Test eta is smooth, spatially Schwartz and compact in I.
Its amplitudes are unrestricted; DS_eta is linear in eta.

The complete source vanishes at the zero scalar jet.
Integrate its first jet derivative along the real segment
from zero to the actual jet. With15 inputs, four outputs
and per-component Cauchy bound M, the Euclidean vector
L2 source norm satisfies ||S||2<=30M U_Phi.
Spatial differentiation gives the same bound for each
spatial derivative. Summing three components by
sqrt(3)<2 gives ||grad S||2<=60M U_Phi.

Likewise ||DS_eta||2<=30M U_eta.
For each spatial derivative of DS_eta, the direct term
has15 first derivatives; the differentiated coefficient
term has15^2 second derivatives, bounded by2M, multiplied
by third Phi jets bounded pointwise by1. Thus each of
its four output components is bounded in L2 by
(15+2*15^2)M U_eta=465M U_eta.
The four-output and three-spatial-component sums give
||grad DS_eta||2<=1860M U_eta.

Set C=30M=600000N/(kappa*kappa0). The four bounds become

    ||S||2<=C U_Phi,       ||grad S||2<=2C U_Phi,
    ||DS_eta||2<=C U_eta,  ||grad DS_eta||2<=62C U_eta.

These are continuous full-function bounds, including the
constant Href construction and every R-dependent term.
The large canonical kappa0 hierarchy, not a deleted
high-frequency tail, supplies their small coefficient.
