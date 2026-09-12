# Full complex exterior coercivity from the unchanged factors

Set zeta=p/(4m^2) and d=1+1/zeta. The source-pinned radial representations define Atrace=-Ftrace and A2=-F2 analytically on the first sheet slit at p<=-4m^2, with removable p0. Their finite coefficients remain4 and1/30. In particular A2 has the same original subthreshold zero and its inverse has the same pole.

For the upper cut p=-tau+i0, put z=1-4m^2/tau. The actual shear real part is

D2=-172/225+19z/30-z^2/10
   +sqrt(z)(30-20z+3z^2)atanh(sqrt(z))/30.

Its polynomial increment is z[8/15+(1-z)/10]>=0, and the last term is nonnegative. Hence D2>=-172/225 on both banks. At the threshold this lower bound is attained. The original S86 trace factor already obeys Re Atrace>=16/15 throughout the first sheet.

## Uniform outer-circle asymptotics, including both sides of the cut

For |zeta|>=2, the principal sqrt(d) lies near1 and has positive real part. The identity

atanh(1/sqrt(d))=[Log(zeta)+2Log(1+sqrt(d))]/2

is fixed by equality on positive zeta and analytic continuation on the slit exterior. It is not a fresh principal-branch reset at each point. The factors sqrt(1+1/zeta) and Log(1+sqrt(1+1/zeta)) are analytic in1/zeta on a disc about0, with uniform Taylor remainders. Substitution into the exact closed forms gives

Atrace=2Log(zeta)+4log2-14/15+O((1+log|zeta|)/|zeta|),
A2=(13/60)Log(zeta)+(13/30)log2-52/225
   +O((1+log|zeta|)/|zeta|),

uniformly on the entire large circle approached from either cut bank. The imaginary part of Log(zeta) is bounded bypi. Production checks the exact leading coefficients; uniformity follows from the analytic convergent coefficient expansions, not from a real-axis asymptotic extrapolation.

Apply the harmonic minimum principle to Re A2 on a large slit disc, excluding a small threshold circle. Both banks and the threshold limit are bounded below by-172/225, while sufficiently large outer circles exceed that value uniformly. Taking the two limits proves

Re A2(p)>=-172/225

on the complete first sheet. Zeros of A2 are not poles of A2 and do not obstruct this harmonic argument.

## A common logarithmic minorant on the slit exterior

Let u=log(tau/(4m^2))>=0 along the cut. For u>=16, z>=15/16 andsqrt(z)>=3/4. The shear polynomial weight is at least13, its polynomial part exceeds-1, and atanh(sqrt(z))=u/2+log(1+sqrt(z))>=u/2. Therefore

D2>=13u/80-1.

For trace, 3-2z+3z^2>=8/3 and the polynomial part16/15+z-3z^2 exceeds-1. Hence Dtrace>=u-1 on the same tail. On0<=u<=16, the global shear and trace gaps give, with exact positive margins,

D_i>=13u/80-4.

Thus the latter bound holds on BOTH complete cut banks, in both channels.

Now apply the minimum principle to
Re A_i(p)-(13/80)log(|p|/(4m^2))
on the slit annulus4m^2<|p|<R. The inner circle is bounded below by-4 using the global gaps; the cut banks have the same bound; outer circles tend uniformly to positive infinity because both actual leading log coefficients,13/60 and2, exceed13/80. Threshold continuity and the limiting inner circle handle their intersection. LettingR grow gives

Re A_i(p)>=(13/80)log(|p|/(4m^2))-4

on the full first-sheet exterior. This continuous complex-plane proof is independent of the numerical diagnostic grid.
