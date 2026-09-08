# Prepared variation at fixed acoustic time

Use physical N=1+e n, a_hat=a exp(e zeta), fixed scalar clock and
the original mass functions. The source norm is the maximum of
both fields and u derivatives zero through two on I=[-1/2,1/2].
They are smooth and zero near u0=-1/2.

At e=0 define

    chi=delta log A=zeta+(alpha+beta)n/4,
    r_sigma=delta log f_L=[1+(beta-alpha)/2]n-zeta,
    j=delta log U=2zeta+alpha n.

Here alpha=4/(9h), beta=28/(81h), h=(1+u^2)^3. The source map
(n,zeta)->(chi,r_sigma) has determinant -(1+4/(27h)), never zero.
In particular n=(chi+r_sigma)/(1+4/(27h)); the inverse introduces
no division by H, the constraint mixing, or spatial momentum.
Its source direction chi equals the one in S6.69's rank-one pole.

Primes in this note denote the original acoustic derivative a d/du.
Varying this derivative gives

    delta(D_sigma F)=D_sigma(delta F)-r_sigma D_sigma F.
    delta U=Uj,
    delta U'= (delta U)'-r_sigma U',
    delta U''=(delta U')'-r_sigma U'',
    delta b=chi'-r_sigma b,
    delta(A''/A)=chi''+2bchi'-b r_sigma'-2r_sigma A''/A.

Direct differentiation of the varied time derivative and pump checks
these identities. The highest source derivatives of delta V_infinity
are -a^2[chi]_uddot; lower source jets and clock derivatives remain.
At the bounce, omitting the r_sigma terms changes the lapse-source
coefficient by 616/81 and the scale-source coefficient by -8.
This is a genuine time-change control, not an optional convention.

## The retarded coordinate shift

The varied acoustic time is fixed by the same initial endpoint.
At fixed physical u,

    xi(u)=delta sigma(u)=integral_{u0}^u r_sigma(v)/a(v) dv.

It vanishes on the initial neighborhood. Comparing potentials at
a fixed value of sigma instead gives

    delta V|_sigma=delta V|_u-xi V',
    delta R_k|_sigma=delta R_k|_u-xi R_k'.

Since a>=1 and 0<1-4/(81h)<1, |r_sigma|<=2||(n,zeta)||_C0.
The physical interval has length one, so |xi|<=2||(n,zeta)||_C0.
This history term is retained, rather than evaluating only the
fixed-u variation. Prepared initial mode data stay unchanged because
all metric jets and both time-coordinate changes coincide on the
initial neighborhood. The old nonzero mixing is not reset to zero.

## Explicit C2 momentum-decaying bound

Put D=k^2+U and A1=bU'+U''/2. Exact differentiation gives

    delta R|_u = V1/D+V2/D^2+V3/D^3,
    V1=delta b U'+b delta U'+delta U''/2,
    V2=-A1 delta U-(3/2)U' delta U',
    V3=(3/2)U'^2 delta U.

Similarly

    R'=H1/D+H2/D^2+H3/D^3,
    H1=A1', H2=-A1 U'-(3/2)U'U'', H3=(3/2)U'^3.

Two independent differentiations verify the rational variation and
history rows. Each Vp/m^(2p) is source linear and contains no source
derivative above second. On the full continuous interval, exact
positive-denominator coefficient majorants give

    ||V1/m^2|| <=3343645513/5308416 * ||source||_C2,
    ||V2/m^4|| <=40047147671/28311552 * ||source||_C2,
    ||V3/m^6|| <=61181640625/50331648 * ||source||_C2,

and

    |H1/m^2|<=18203125/16384,
    |H2/m^4|<=830078125/262144,
    |H3/m^6|<=91552734375/33554432.

Because U>=m^2, theta=m^2/(k^2+m^2)<=1 and
m^(2p)/D^p<=theta^p<=theta. Thus the fixed-acoustic-time
variation satisfies

    ||delta R_k||_C0 <=
      (11735921060393/679477248) theta ||source||_C2
       <18000 theta ||source||_C2.

The factor two multiplying the three history bounds comes from
the actual retarded xi estimate. No massless limit or k=0
canonical coordinate is used to establish it. The potential-level
bound extends to k=0 by continuity; the underlying canonical
transformation remains stated for k>0.
