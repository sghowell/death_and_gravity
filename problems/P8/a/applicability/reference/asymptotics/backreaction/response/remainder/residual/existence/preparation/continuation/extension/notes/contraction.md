# Complete weighted contraction with the original short preparation

All comparisons are on t=x-x0 in [0,L], L=10^-6. The source still uses
L0=10^-10. Let sigma=2*10^6, alpha=1/sigma, E=9, r=2*10^-7 and
`||v||_sigma=sup exp(-sigma t)|v(t)|`. The ball is X(0)=0,
||X||_sigma<=r, with u=ubar+I X. The past and quantum modes are common.

## 1. Input provenance and pointwise enclosure

The unchanged A.8/A.11 calibration supplies

    |ubar|<=delta b0, b0=20461/128,
    |ubar'|<=delta b1<V=10^-10, b1=61013499/8192,
    |qbar|,|qbar'|<=Qb=10^-8,
    |cbar|<C=12*10^-9.

These bounds hold on the old plateau, not just its previous short slab.
Its exact clock has dy/dx<=1, so y0+L<3 keeps the baseline in that
plateau. The actual original-history duration is still <=T=3. The
original source is c=chi_L0*cbar; consequently |c|,|c-cbar|<=C on the
whole longer interval, although c is already zero from L0/2 onward.
No derivative bound on the cutoff or cbar is needed in C0.

For any weighted ball element, causality gives

    |X(t)|<=E r, |I X(t)|<=E r alpha,
    |u|<=delta b0+E r alpha < B=3*10^-12,
    |u'|<=V+E r < M=10^-5.

The unchanged past has the same M bound. The pointwise exponential is
not ignored: sigma L=2 and e<3 prove exp(2)<E=9. For an independent
rational witness, the exponential series through degree three plus its
geometrically bounded tail gives e<87/32<3.

The old exact initial values obey 12/5<a0<13/5 and 3/8<h0<5/12.
While 1/3<=h<=1/2, h'=-u-h^2 obeys |h'|<1/3. The first-exit margins
L/3<1/24 keep h in that box. Since a'=ha and h>0,

    2<a<= (13/5)exp(L/2)<(13/5)/(1-L/2)<3.

These margins, including the exact fourth-power initial inequalities,
are recomputed for the new length, not copied from the old time bound.

## 2. Weighted metric gains and the local frozen coefficient

For two inputs write D=||X1-X2||_sigma. The primitive has norm <=alpha.
The same-data Hubble difference satisfies the exact damped equation

    Delta h(t)=-integral_0^t exp[-integral_s^t(h1+h2)] Delta u(s)ds.

The damping is <=1. Applying weighted integration repeatedly and then
the mean-value theorem on 2<=a<=3 gives the constants, per D,

    Delta u: alpha, Delta h: alpha^2, Delta log a: alpha^3,
    Delta a^2:18alpha^3, Delta a^-2:alpha^3/2,
    Delta d:alpha^3/2, Delta d':alpha^2/2.

All are weighted norms. For the fixed af=5/2,
`a0^4=af^4-delta`, so |a0^2-af^2|<=delta. Also |(a^2)'|<=9.
Thus `|a^2-af^2|<=Ea=delta+9L`. Since

    |d(a0)-d(af)|<=delta/[8(af^4-delta)]<delta,
    d'=-h/2,

we have `|d(a)-d(af)|<=Df=delta+L/4`. This small coefficient is
available only after moving the full d(af) term into the A.13 inverse.

## 3. Full auxiliary equation and source terms

Use precisely the same initial q,q' and full auxiliary equation as A.11:

    P=a^2q',
    P=P0+I[a^2u/(60delta)-Pcurv(h,u)]
       +8[c/h-c(0)/h0-I(c(1+u/h^2))],
    q=q0+I(P/a^2),
    Pcurv=u^2/4+h^2u/30+h^4/30.

Both source endpoints are included. The identity `(1/h)'=1+u/h^2`
removes c' exactly; using a radiation pressure while switching c would
not be the same conserved source.

At X=0, the geometry and actual Wick function are the old ones. Its
source is cbar, whereas the new auxiliary source is c. Therefore

    EP=8C[3+L(1+9B)],
    ||Pcenter-Pbar||<=EP,
    ||qcenter-qbar||<=L EP/4,
    ||Gf[0]||_sigma<=EG=EP(1+9L/4).

The source is not stretched to L in this estimate: only |c-cbar|<=C
is used. All other actual-state/local differences vanish at the center.

For a pair of candidates, weighted mean-value bounds give

    KC=(B/2+1/120)alpha^2+(B/30+1/60)alpha^3,
    Ksource=144C alpha^2+432CB alpha^3,
    KEfull=(9alpha^2+18Balpha^4)/(60delta),
    KP=KEfull+KC+Ksource.

For the source, `Delta(1/h)<=9Delta h` gives the 72C endpoint term;
`Delta(u/h^2)<=9Delta u+54BDelta h` gives the other 72C primitive
term and 432CB term. The source function itself is fixed in the pair.

Before using any proposed auxiliary caps, verify without circularity

    |P|<=9Qb+EP+E r KP <Pmax=10^-5,
    |q|<=Qb+L Pmax/4 <Qmax=2*10^-8.

KP is independent of Pmax,Qmax. Its leading value is 3.75, but it
multiplies the small radius; the derived P cap is below 7.129*10^-6.
Now weighted integration gives

    Kq=alpha KP/4+Pmax alpha^4/2,
    Kaux=9Kq+18Qmax(alpha^2+alpha^3),

where Kaux bounds the complete difference of 2a^2hq.

## 4. The full resummed RHS

Subtracting the baseline, A.13's exact decomposition becomes

    (Lf-cf I^2/2)X=Gf[X], cf=af^2/(30delta),
    Gf=rbar+I[(a^2-af^2)u-(abar^2-af^2)ubar]/(60delta)
         -I Delta Pcurv+Delta(source primitive)
         +2Delta(a^2hq)-Delta N',
    N=R[u]+(d(a)-d(af))u.

Here rbar is the actual center defect, not set to zero. No quantum
history, anomaly, auxiliary or source term is deleted. The remaining
Einstein coefficient has pair constant

    KErem=(Ea alpha^2+18Balpha^4)/(60delta).

For the non-mode part of N', use the exact pair decomposition

    Delta[(d-df)u']=(d1-df)Delta X+Delta d*u2',
    Delta[d'u]=d1'Delta u+Delta d'*u2.

Its full bound is

    Klocal=Df+(V+E r)alpha^3/2+alpha/4+Balpha^2/2.

The actual nonlinear response is bounded by the A.10 shared-history
lemma. At each observation t it gives a coefficient

    KR(t)=MT t^2[5/4+log(T/t)/2]+18M^2T^5 t exp(2MT^3)

times sup_[0,t]|Delta X|. The quadratic coefficient extends continuously
by zero to t=0; its derivative is `t[2+log(T/t)]>0` for 0<t<=L<=T.
The other coefficient is increasing too. Since
`sup_[0,t]|Delta X|<=exp(sigma t)D`, the weighted pair bound is KR(L),
**without** an extra exp(sigma L) factor. A rational coarsening is

    KR<=MT L^2(5/4+10)+36M^2T^5L.

Here log(T/L)<20 follows, for example, from (8/3)^20>T/L and e>8/3;
exp(2MT^3)<2 follows from 2MT^3<1/2. The code checks these assumptions.
The modes retain the full original shared history; there is no state reset.

The sum of every displayed contribution is

    KG=KErem+KC+Ksource+Kaux+Klocal+KR
      =25031252674164079302312578725005874194061
        /3000000000000000000000000000000000000000000000
      <9*10^-6.

## 5. Actual Banach gate and pointwise conversion

A.13 proves the complete causal inverse has weighted norm
`k<=240/769` at sigma=2*10^6, including its positive pole. Hence

    q=k KG<3*10^-6,
    eta=k EG<9*10^-8,
    ||C[X]||_sigma<=eta+q r<r=2*10^-7,
    ||X*||_sigma<=eta/(1-q)<9*10^-8.

The last strict margin uses the exact eta and q, not just their coarse
individual upper bounds above. Every inequality has a strict rational margin. This is the complete
map's self-map/contraction, not a conditional number assigned to an
unbounded omitted term. Converting back with E=9 proves
`||X*||_infinity<81*10^-8`; every iterate has the earlier pointwise
metric/state/auxiliary barriers. Smoothness, the actual stress tensor,
initial constraint, and agreement with A.11 follow in the companion proof.

The new duration is not macroscopic continuation and does not suppress
the growing pole. No physical prescription or source support was changed
to obtain the bound. It uses the still modest exp(2) conversion, which
cannot be assumed harmless on arbitrary larger intervals.
