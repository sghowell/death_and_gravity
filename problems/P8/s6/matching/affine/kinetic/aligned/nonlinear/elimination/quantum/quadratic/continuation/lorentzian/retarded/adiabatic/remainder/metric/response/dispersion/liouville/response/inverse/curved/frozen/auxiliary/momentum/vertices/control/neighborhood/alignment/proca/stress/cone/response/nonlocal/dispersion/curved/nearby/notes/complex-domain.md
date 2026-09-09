# A fresh joint complex-time lapse root

The upstream S6.75 Hamiltonian majorant itself was proved for real u,
though its coefficients and boundary primitive were analyzed on a
complex-time tube. A bound on I_u at real u cannot simply be assigned
to every complex u in that same tube. This new proof uses a smaller
time disc and recomputes the needed derivative bound.

## Complex coefficient domain

For each central N0 in [1,1+10^-6], use

    |u|<=1/100, |N-N0|<=1/200,
    |p|<=1/4, |ell-ell0|<=1/4.

The lapse disc lies strictly inside the original |N-1|<=1/100
domain. The upstream literal primitive is holomorphic and bounded
by 2/99 on |u|<=1/50, for the same complex lapse domain.
A circle of radius 1/100 about any u in the smaller disc remains
in that primitive domain. Cauchy's formula therefore gives the
fresh joint bound

    |I_u| <=200/99<3.

Strict nonzero denominator margins give an open neighborhood of
these closed domains, so no unjustified boundary differentiation
is needed. Principal powers use their original common branch.

The original coefficient estimates for this smaller domain give
|U|,|U^-1|<2, |(4a)^-1|<2 and |b|<1. The full transformed scalar
plus margin remains below 1002. Consequently

    |f0| <2*1002+2*(200/99)<3000.

Since ell0<11/100, |ell|<9/25 in the outer disc. The homogeneous
Hamiltonian is therefore bounded by

    |Hcal| <2[2*(5/4)^2+2*1002+2*(200/99)+(9/25)^2]
           <10000.

This is a bound on the actual full homogeneous Hamiltonian with
its primitive, not on a finite central Taylor polynomial. The
source-free mass update is exact and the zero-vector background
does not contribute any hidden source term.

## Recentered contraction

On the inner half-polydisc, remaining Cauchy radii are
r_N=1/400, r_u=1/200 and r_p=r_ell=1/8. With M=10000,

    |H_Nu|<=800000000,
    |H_Np|,|H_Nell|<=32000000,
    |H_NNN|<=3840000000000,
    |H_NNu|<=640000000000,
    |H_NNp|,|H_NNell|<=25600000000.

The factors 6 and 2 in the third-derivative bounds are retained.
At the exact central fixed phase, F=H_N=0 and h0=H_NN<-2,
uniformly by acceleration.md.

Take r=10^-24 for u,p,ell-ell0 and RN=10^-14 for N-N0.
The straight-segment estimate gives

    |F_N-h0|<=3840000000000*RN
              +(640000000000+2*25600000000)*r=:d,
    |F(u,p,ell,N0)|<=(800000000+2*32000000)*r=:f.

Use the fixed-center map n -> n-h0^-1 F(u,p,ell,N0+n).
Its Lipschitz constant is at most d/2<1/10, and its center
step is at most f/2. The exact check has

    f/2+(d/2)*RN<RN.

Thus the closed complex lapse disc is invariant and the map
is a strict contraction, uniformly on the entire parameter box.
The root is unique in this disc and satisfies

    |N-N0|<=f/(2-d)<10^-15,
    |H_NN|>=2-d>19/10.

Uniformly convergent holomorphic iteration gives a holomorphic
root in u,p,ell; strict margins extend it to an open neighborhood
of the parameter box. Conjugation symmetry and uniqueness make
it real for real inputs. Its lapse is positive and its X=N^-2
stays in the original covariant clock tube.

This is newly centered at the chosen finite matter-charge shift.
The larger separation ell0-1/10 is not assigned the old
10^-24 nine-invariant polydisc centered on the original clock.
The old action and old certified domains are unchanged.
