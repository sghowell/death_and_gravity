# Actual time derivatives and differentiated subtraction

For each actual transverse/longitudinal energy or pressure weight
A_Q,B_Q, write p=v'-d v and

    Q=(A_Q |p|^2+B_Q omega^2 |v|^2)/(2a^3).

The exact mode equations give v'=p+d v, p'=-omega^2 v-d p.
Consequently X=(|v|^2,Re(v conj(p)),|p|^2)^T obeys X'=M X,
where

    M=[[2d,2,0],[-omega^2,0,1],[0,-2omega^2,-2d]].

This identity is independently checked using the four real canonical
coordinates. With r_0=(B_Q omega^2,0,A_Q), define

    r_{j+1}=D r_j-3H r_j+r_j M.

Then Q^(j)=r_j X/(2a^3). The full fixed-comoving derivative D
acts on u,z and omega^2; no truncated Hubble-jet rule is used.
For j=0,...,5 every coefficient is a polynomial in omega^2
with rational u,z coefficients. Exact boxes give E_j such that
sum |r_j,c|*(omega^-1,2,3omega)_c <= E_j omega^(j+1).
All powers satisfy this inequality because omega>=m>=1000>1.

## Compare products, not successive approximate equations

The exact CCR and |Bcal|<1 give changes in each reference
squared product at most 6|Bcal| times that square, and in the
mixed real product at most 6|Bcal| |f p_f|. Thus

    |Q^(j)(v)-Q_projected,j(f)|
        <=3 E_j |Bcal| omega^(j+1)/a^3.

Here Q_projected,j means the exact-ODE row r_j evaluated on
the reference products. It is generally NOT the j-th derivative
of Q(f): f has a residual. Defining the projection explicitly
avoids dropping those residual terms while differentiating.

For W=omega*S, R=S', c=c1+R/(2S), the projected reference
equals omega*F_j/(4a^3), with

    F_j=[t(r_X-c r_Y+c^2 r_Z)+S^2 r_Z]/S.

The frozen fourth-order subtraction bracket has coefficients
F0,F2,F4. Its actual j-th derivative is obtained coefficientwise
by j applications of D0+[(1-2n)lambda-3H] to coefficient n.
Call this bracket AD_j. An independently checked denominator
identity is

    4t^2 S^3 (F_j-AD_j)
      =t^2{4t(r_X-c1 r_Y+c1^2 r_Z)S^2
          +2t(-r_Y+2c1 r_Z)R S+t r_Z R^2
          +4r_Z S^4-4AD_j S^3}.

The right side is an ordinary polynomial in t for j<=5. For
the eighth-order reference, its coefficients of t^0,...,t^4
vanish exactly for all 24 actual readout/order pairs.
Replacing every remaining coefficient operation by its absolute
polynomial envelope, and using 1/(4S^3)<=2, gives

    |F_j-AD_j|<=T_j t^3,
    |Q_projected,j-AD_density^(j)|<=T_j/(4a^3 omega^5).

The report gives the individual integer T_j; their largest is
4275264710183533. This is a conservative continuous envelope,
not a sampled maximum.

A failure control repeats the calculation with only W_4 at
j=2 for longitudinal energy. The t^4 numerator coefficient
at u=0,z=1 is 31744/81, not zero. Thus that lower reference
does not support this separate absolute-tail argument. This
does not establish divergence of the exact state's derivative:
separately nonintegrable comparison terms could cancel.

## Continuous integrals and physical normalization

Let K=max(8m,10^12). After summing three polarizations, use
omega<=A*nu and nu^(j-5)<=m^(j-5). The radial measure is
k^2 dk/(2pi^2). At the maximal power j=5, the initial bands
and all-momentum evolution give

    integral |Bcal| nu^6 d^3k/(2pi)^3 <= J5,
    J5=A^3[32 B6 m^3/27+B8 K/18+B10/(9K)+Kmix/(24m)].

For the low band, k<=4Am and pi^2>9 bound its volume.
For the middle band k^2/nu^2<=A^2 and k<=AK.
For the high band k>=A sqrt(K^2-m^2)>=AK/2 and
k^2/nu^4<=A^4/k^2. The last contribution uses the exact
integral of nu^-4, A^3/(8pi m), followed by pi>3.
Its radial antiderivative is independently differentiated.
All these bounds hold for the full continuum, not a grid.

It follows that the differentiated state/reference difference
for either observable is bounded by 9 E_j A^(j+1)m^(j-5)J5.
The projected subtraction tail integrates to T_j/(72m^2),
because the a^3 from integrating omega^-5 cancels a^-3.
Differentiate the actual S6.53 matched local coefficients
directly j times in u and bound them by L1 polynomial boxes.
Their contribution is sum_n U_{n,j}m^(4-2n)/(576L^2),
where L=M*tau and the denominator uses pi^2>9.

Divide the two state contributions by L^2 as well. The result
bounds the j-th physical-time derivative divided by
M^2/tau^(2+j). At L=10^24,m=1000 every bound for j<=5 is
<10^-18; the fifth-order energy and pressure bounds are,
respectively, <1.165*10^-19 and <1.132*10^-19.
The certificate stores exact rationals, not these rounded displays.
The example changes a scale parameter only: it does not overwrite
the earlier L=10^12 example or alter the state or counterterms.

The uniform integrable envelopes for each actual derivative of
the exact-minus-subtraction integrand justify differentiation
under the integral through order five, with continuous one-sided
endpoint limits. Hence the finite matched stress is quantitatively
C5 on I. This is fixed-background regularity, not a bound on
functional metric variations, compatible-state changes, a full
renormalized causal feedback map or a quantum-corrected bounce.
