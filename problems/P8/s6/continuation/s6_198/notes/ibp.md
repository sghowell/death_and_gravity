# Exact retarded identity, including every endpoint

For one full reference polarization pair and fixed k,l, write

    K(t)=exp(-iTheta(t))*integral_(t0)^t exp(iTheta(s))b(s)ds,
    Theta'=Omega=W_k+W_l, g=1/Omega, Lb=(g b)'.

Here b is the complete phase-stripped reference tensor
amplitude contracted with Gammahat(t,k+l).
It is not a scalar replacement for the physical vertex.

Since (exp(iTheta))'=iOmega exp(iTheta), repeated integration
by parts gives, for n<=6,

    K=B_n+i^n exp(-iTheta(t))*integral_(t0)^t exp(iTheta(s))L^n b(s)ds,
    B_n=sum_(j=0)^(n-1) (-i)i^j g L^j b.

Every lower endpoint is zero: Gamma and all its jets vanish
in a neighborhood of the unchanged preparation surface.
If that hypothesis is removed, lower endpoints must be kept.

The direct first-order equation is K'+iOmega K=b.
The exact algebra checks

    B_n'+iOmega B_n=b-i^n L^n b

for all six orders. The remainder supplies precisely the
missing term and has zero common initial value.
The six boundary coefficients are

    -i, 1, i, -1, -i, 1.

The reference current is minus the imaginary part after
contracting with the complete conjugate detector pair and
summing all9 polarization pairs. Complex amplitude terms
can make odd endpoints nonzero. No such term is removed
before a full physical tensor contraction demonstrates it.

In particular B_6-B_5=g L^5 b. This final endpoint has six
inverse-phase factors and is separately integrable together
with the sixth bulk. The remaining first five endpoints
are left explicit for fixed spatial matching.
