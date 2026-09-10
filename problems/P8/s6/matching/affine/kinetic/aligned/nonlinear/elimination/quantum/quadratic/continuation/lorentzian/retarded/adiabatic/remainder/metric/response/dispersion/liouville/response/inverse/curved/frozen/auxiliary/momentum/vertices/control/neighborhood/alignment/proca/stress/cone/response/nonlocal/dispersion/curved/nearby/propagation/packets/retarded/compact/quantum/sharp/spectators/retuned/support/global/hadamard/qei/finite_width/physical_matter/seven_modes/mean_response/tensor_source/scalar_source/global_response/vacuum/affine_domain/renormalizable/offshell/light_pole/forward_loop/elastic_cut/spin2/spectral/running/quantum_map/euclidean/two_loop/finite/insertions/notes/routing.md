# A common first-sheet routing for all three channels

In bilinear Euclidean coordinates choose

    p1=(i sqrt(4-s)/2,0,0,i sqrt(s)/2),
    p2=(-i sqrt(4-s)/2,0,0,i sqrt(s)/2),
    p3=-p1, p4=-p2.

The square roots are analytic on |s-2|<=1. Each pi^2=-1;
the three pair invariants are s, 0 and 4-s. Every identity
is checked literally, including momentum conservation.

Hermitian norms are used only for estimates, not to change
the bilinear momentum invariants. On this disc,

    |pi|^2=(|s|+|4-s|)/4<=3/2,
    |P_channel|^2<=3.

In an outer bubble the internal lines are k and k+P.
Feynman combination and the shift q=k+xP give

    [q^2+Delta]^-2,  Delta=1-x(1-x)s_channel.

Its real gap is at least 1/4. The inserted self-energy
lines are q-xP and q+(1-x)P; their shifts have squared
Hermitian norm at most three, as used in notes/kernel.md.

At the first full quartic vertex the two internal heavy
exchange shifts are pa-xP and pb-xP. Their norms are at
most sqrt(3/2)+sqrt(3)<3. The opposite vertex obeys the
same bound by momentum conservation. After a real
translation, each full heavy denominator has real part
at least M-9>M/2. The external-pair heavy exchange has
the still stronger lower bound M-3.

Hence both full vertices satisfy

    |V| <= lambda4+6g/M < 3lambda4

at the actual parameters. This retains the local term
and all three heavy channels for each vertex.

The identities first hold with real Euclidean external
momenta and a common translation-invariant regulator.
Scale the final complex momenta from zero to obtain the
specified first sheet. All displayed margins persist
along this path: the combined light denominator,
heavy denominators and auxiliary kernel denominators
remain separated from zero. The decaying-remainder
integral is absolutely dominated, so its continuation
and contour shift survive removal of the regulator.
The constant piece instead uses the already subtracted
S6.113 amplitude and is never treated as a finite bare
integral.

The strict margins also permit a neighborhood of the
closed unit disc. Cauchy's estimate can therefore be
applied to the integrated remainder, not only formally
to its pointwise integrand.
