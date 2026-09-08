# Full Proca two-insertion pole

## Fixed component and conventions

Keep S6.42 unchanged, with canonical vector Wc=sqrt(zeta_physical)W
and m2=m0^2=M^2/zeta_physical. Freeze the physical metric and the
on-clock mass tensor at a=b=1. The first mass variation is
delta M=m2 A r(x), A=diag(alpha,beta,beta,beta),
alpha=4/(9h), beta=28/(81h), h=(1+u^2)^3. These numbers are
derived directly from the original gamma_t,gamma_s and p_N=-1/(2h).
The label r varies this tensor with the metric held fixed. It is not
the full physical lapse fluctuation of the gravitational theory.

Use a Euclidean +1/2 Tr log K and d=4-2epsilon_DR. For
p^2=sum p_mu^2 the full propagator is

    G(p)=[I+p p^T/m2]/(p^2+m2),
    K(p)=(p^2+m2)I-p p^T.

Their product is the identity. In particular, the longitudinal
numerator is not discarded as a supposed gauge mode. For a general
real symmetric A the two-insertion trace numerator is exactly

    m2^2 tr(A^2) + m2[p^T A^2 p+q^T A^2 q] + (p^T A q)^2.

The checker expands all ten independent entries of A. Expanding
the Gaussian logarithm gives -1/4 Tr(G delta M G delta M). Therefore,
if F denotes the integrated numerator's pole weight, the quadratic
functional is

    Gamma_bubble,pole^(2) = -1/(64*pi^2*epsilon_DR)
       integral d^4k/(2*pi)^4 r(-k)r(k) F(k).

There is no implicit extra one half in this convention. If instead
Gamma^(2)=1/2 integral r(-k)Pi(k)r(k), Pi's pole is -F/(32*pi^2*epsilon_DR).

## Angular and radial derivation

Combine the two denominators with x in [0,1], t=x(1-x), q=p+k,
p=l-xk and q=l+(1-x)k. Write Delta=m2+t K,
K=k^2, T=tr A, U=tr(A^2), B=k^T A k and C=k^T A^2 k.
The pole-level angular moments in four dimensions are
<l_i l_j>=l^2 delta_ij/4 and
<l_i l_j l_k l_l>=l^4(delta_ij delta_kl+delta_ik delta_jl+
delta_il delta_jk)/24. Relative to 1/(16*pi^2*epsilon_DR), the
radial residues for 1,l^2,l^4 over (l^2+Delta)^2 are
1,-2Delta,3Delta^2. Actual Gamma residues are checked separately.
Dimensionally continued angular coefficients affect finite local
terms, but not the pole evaluated at d=4.

The resulting parameter integrand is

    F(x,k)=m2^2 U-m2 Delta U+m2[x^2+(1-x)^2]C
           +Delta^2(T^2+2U)/8-(1-2x)^2 Delta C/2
           +t Delta T B+t^2 B^2.

Its exact x integral splits into

    F0=m2^2(T^2+2U)/8,
    F2=m2[K(T^2-2U)/24+C/2+T B/6],
    F4=K^2(T^2+2U)/240-K C/60+K T B/30+B^2/30.

An independent implementation expands the shifted four-component
polynomial at k=(k_time,k_space,0,0), applies explicit spherical
monomial averages and then the radial residues. For the rotationally
symmetric mass direction this covers arbitrary real momentum.
Polynomial identities then also continue to complex momentum.

## Actual coefficient and continuous nonvanishing

For T0=k_0^2 and S0=sum_i k_i^2 at real Euclidean momentum,

    F0=904*m2^2/(2187*h^2),
    F2=416*m2*(4*T0+3*S0)/(6561*h^2),
    F4=32*(120*T0^2+220*T0*S0+101*S0^2)/(98415*h^2).

Since 2*101<220<2*120, its numerator lies between
101*(T0+S0)^2 and 120*(T0+S0)^2. On |u|<=1/2,
1<=h<=125/64, and the exact endpoint comparisons give

    (k^2)^2/125 < F4 < (k^2)^2/25

for k nonzero. This is an actual, nonzero fourth-derivative pole of
this off-shell mass component. It is not the pole of a constrained
physical scalar propagator or evidence of a new physical ghost.

## Independent reductions and boundaries

For A=c I, Dp=p^2+m2, Dq=q^2+m2, use
p.q=(Dp+Dq-2m2-K)/2. Direct rational algebra gives

    N/(Dp Dq)=(3m2^2+m2*K+K^2/4)/(Dp Dq)
                -K(1/Dp+1/Dq)/2+(Dp/Dq+Dq/Dp+2)/4.

A loop shift reduces integral(Dp/Dq) to K times a tadpole, up to
odd and scaleless pieces. The scalar bubble residue is 1 and the
tadpole residue is -m2. This independently gives
F=c^2(3m2^2+3m2*K/2+K^2/4), in agreement above.

At zero momentum, the coefficient of r^2 in S6.47's pole weight
with a=1+alpha*r and b=1+beta*r is
[(alpha+3beta)^2+2(alpha^2+3beta^2)]/8=F0/m2^2.
This uses a linear mass variation. The original mass tensor also
has a nonzero second variation whose tadpole contributes to the full
clock second derivative. A negative control rejects identifying the
bubble alone with that full derivative. Another rejects omitting the
longitudinal numerator, which loses all momentum-dependent poles.

These derivations fix the regulator and Gaussian sign directly; no
sign-normalized agreement with a differently printed epsilon convention
is assumed. Counterterms, full metric/scalar variations and other field
loops are not determined by this component. The finite result and its
strictly limited domain are proved in [the remainder note](remainder.md).
