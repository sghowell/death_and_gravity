# Whole dimensional MS reference and its exact finite conversion

Use d=4-2epsilon, tr 1=4 and exp(gamma_E epsilon) mu^(2epsilon)
per loop, exactly as the inherited scheme. Set CY=2NY^2/Q.
Before setting mu=m, the second mass derivative of the regulated
zero-momentum two-point kernel gives

 Gamma0_D(0)/(CY)=(mu^2/m^2)^epsilon
 exp(gamma_E epsilon) Gamma(epsilon)(12-32epsilon+16epsilon^2).

The proper box counterterm is -12CY/epsilon; its MS finite anchor is
-32CY. Do not apply the MS operation only after multiplying a
four-dimensional finite vertex by the divergent outer light bubble.

For 0<Re epsilon<1/2 the raw zero-soft kernel has the convergent
unsubtracted spectral representation -integral nu_D(v)/(v+q^2) dv.
It follows by differentiating the regulated two-point representation
at fixed mu, or by its Feynman parameters and vanishing large-q
limit in this strip. With p=3/2-epsilon,

 nu_D/(CY)=exp(gamma_E epsilon) mu^(2epsilon) 4^epsilon sqrt(pi)
 /[2Gamma(3/2-epsilon)] *8p v^(-epsilon)
 *(1-T/v)^(-1/2-epsilon)*[(2-2epsilon)T/v-1].

The endpoint powers are integrable in the stated strip. The complete
proper-subgraph-paired outer reference is

 F_D=-integral_T^infinity nu_D(v) T3_D(v,1,1) dv
      -(12CY/epsilon) I2_D(1),

 T3_D=exp(gamma_E epsilon) mu^(2epsilon) Gamma(epsilon)/Q
      *[1/(v-1)-(v^(1-epsilon)-1)/((1-epsilon)(v-1)^2)],
 I2_D(1)=exp(gamma_E epsilon) mu^(2epsilon) Gamma(epsilon)/Q.

For the leading large-v reference replace the bracket by
1/v-v^(-1-epsilon)/[(1-epsilon)v]. Keep both terms: discarding the
first by setting the light mass to zero would erase the soft region.
With z=T/v, the two spectral beta moments are

 I_n=B(n epsilon,1/2-epsilon)
     *[(2-2epsilon)n epsilon/((n-1)epsilon+1/2)-1], n=1,2.

At mu=m and ell=log(m^2), the soft region plus proper counterterm,
in units CY/Q, is

 exp(ell epsilon)[exp(2gamma_E epsilon)(12-32epsilon+16epsilon^2)
 Gamma(epsilon)^2 -12 exp(gamma_E epsilon)Gamma(epsilon)/epsilon].

Its poles/finite term are -32/epsilon+16+pi^2-32ell.
The hard piece is

 -8sqrt(pi) exp(2gamma_E epsilon)4^(-epsilon)
 (3/2-epsilon)(1-4epsilon)Gamma(epsilon)Gamma(2epsilon)
 /[(1-epsilon)(1+2epsilon)Gamma(1/2+epsilon)].

It equals -6B(epsilon)/epsilon^2. The exact logarithmic jets of B
are B(0)=1, (log B)'=-17/3, (log B)''=-103/9+pi^2/3.
Its Laurent terms are -6/epsilon^2+34/epsilon-62-pi^2.
Thus the combined leading reference has poles -6/epsilon^2+2/epsilon
and finite coefficient -46-32ell. The pi^2 cancellation is derived,
not imposed; the finite pole products on both sides are retained.

The difference from this leading reference is uniformly analytic near
epsilon zero (for example |epsilon|<1/8): at large v it gains a full
inverse power and at threshold retains an integrable power. The
removable Gamma(epsilon) singularity can therefore be taken inside.
At d=4, Q*T3=(v log(v)-v+1)/(v-1)^2. Its difference from
(log(v)-1)/v is

 (1/v)[h(k)(-log k-1)+k A(k)], k=1/v,
 A(k)=(1-k)^(-2), h=A-1.

Consequently the exact remaining finite reference is

 delta F/(CY/Q)=-12 integral_0^1 (2z-1)/(z sqrt(1-z))
      *[h(rz)(-log(rz)-1)+rz A(rz)] dz, r=1/T.

For r<=1/16, h(k)<=3k and A(k)<2. Bound the square bracket in
absolute value by k[3(-log k)+5], use |2z-1|<=1,
integral dz/sqrt(1-z)=2, and
integral -log(z)/sqrt(1-z) dz=4(1-log 2)<4.
This gives r[264+72log T], conservatively below r[300+100log T].
There are no remaining pole products in this already convergent
difference. Its nonzero finite contribution is bounded, not set to zero.
