# An exact-frame representative, not an exact quantum solution

For the computed leading mean functions define

    N_eta=1+n,
    hat_a_eta=a exp(xi),
    a_eta=hat_a_eta [(h-1+N_eta^-2)/h]^(-1/4).

This uses the ACTUAL spatial frame map, not its linear
exponential approximation. The linear physical
log-scale coefficient is still B=xi+n/(2h).
For eta<=10^-20 the global bounds give |n|<10^-10.
Since h>=1, the logarithm in the frame factor obeys
|omega|<=|log N_eta|/2<=|n|. Hence
|log(a_eta/a)|<=|xi|+|n|<3*10^9 eta<10^-10.

In particular N_eta>1/2 and a_eta>1/2 globally.
The smooth state source and regular forced ODE give
smooth n,xi at every finite u. Conserved spatial
geodesic momentum then gives divergent timelike
proper length on each tail:

    integral N_eta du/sqrt(1+P^2/a_eta^2)
        >= integral du/[2 sqrt(1+4P^2)].

Null affine length is proportional to
integral N_eta a_eta du, also divergent.
Thus this explicitly defined metric is complete.

Parity makes its physical Hubble parameter vanish at
u=0. The exact spatial-frame second logarithmic
derivative there is

    (log a_eta)''
       =4+xi''-(3/2)[(1+n)^2-1]+n''/[2(1+n)].

The independent band jets give xi''>=0 and n''>=0,
while 0<=n<10^-10 at the anchor. The numerator is
therefore above 4-3*10^-10-(3/2)*10^-20>3.
Dividing by N_eta^2 leaves strictly positive proper
bounce acceleration. This argument keeps the full
frame factor rather than silently truncating it.

These are geometric properties of a representative
of the computed FIRST-ORDER profile. No quantum state
has been co-evolved on this representative, no common
absolute vacuum source has been removed, and no
nonlinear residual bound or exact SEE shadowing has
been established. Therefore this is not an additional
self-consistent quantum bounce or a UV witness.
