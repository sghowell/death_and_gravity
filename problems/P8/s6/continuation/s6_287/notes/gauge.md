# Complete zero-transfer gauge-fixing insertion

Write k=p-r, D_r=r^2-mu and p^2=mu only after applying the literal scalar
Ward identity. The harmonic graviton projector gives
H=P_D T=T-eta*tr(T)/(D-2), with
tr(H)=-2tr(T)/(D-2).
Consequently, for gamma=H/k^2,
F_nu=k^mu gamma_mn-k_nu tr(gamma)/2
     =k.T/k^2
     =[(p^2-mu)r-D_r p]/k^2
     =-D_r p/k^2 on shell.
The trace cancellation is algebraic in D, not a dimension interpolation.

Hold covariant Fourier momenta and the internal covariant gamma fixed
when varying a constant symmetric metric B. Then
delta(g_inverse)=-eta B eta,
delta F=H delta(g_inverse) k_cov
        -k_cov tr(delta(g_inverse)H)/2,
with the common k^-2 supplied by gamma. The whole density variation is

delta(sqrt(-g) F.g_inverse.F)
 =tr(eta B) F.eta.F/2
  +F.delta(g_inverse).F+2F.eta.deltaF.

Every term has at least one F; this is why retaining only a TT
variation or omitting deltaF would not be an adequate argument.
Substitute the proven on-shell F and divide by the internal scalar
propagator denominator D_r. After suppressing the common k^-4, the
complete quotient is the polynomial

D_r[mu*tr(eta B)/2-p.B.p]-2p.deltaF_numerator.

The literal T(p,p-k) is constant plus linear in k. D_r=k^2-2p.k.
The displayed quotient therefore has soft degrees1 and2 and no degree0.
Its odd part integrates to zero. Its even part is homogeneous degree2,
so the remaining loop is a massless quadratic tadpole
integral polynomial_degree2(k)/(k^2)^2. It has no mass or external
denominator and vanishes in dimensional regularization.

Importantly, its radial origin is proportional to
integral_0 dk k^(D-3), convergent at D near4. There is no hidden
degree0/k^4 logarithmic infrared integral whose zero would disguise a
UV/IR pole cancellation. Local quadratic scaleless tadpoles introduce
no dimensional mass scale or finite charge counterterm here.

The code verifies every symmetric metric direction for generic
component momenta in D4,D5,D6, including the harmonic identity, direct
metric derivative, scalar-denominator cancellation, absence of a
constant soft term and degree2 of the whole even numerator. The
general-D proof is the preceding trace and degree argument; the finite
component checks verify the implementation, not an all-D extrapolation.

It follows that the specified ordinary on-shell proper vertex equals
the background one at q0 at this loop order. The background variation
of sqrt(-g)K(g^mn p_m p_n) yields
Gamma0=2pp Kprime-eta K. This equality at the point is not by itself
continuity of the ordinary form factor. That separate obligation is
proved in continuity.md. No off-shell or arbitrary-gauge equality is
claimed.
