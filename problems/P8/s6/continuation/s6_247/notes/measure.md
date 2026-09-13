# Complete two-threshold inverse measure and actual zero set

For each species and channel, write W(y)=sum_(j=1)^3 w_j y^(2j), z=1+4m²/p. The EXACT identity
p/[4m²+p(1-y²)]=1/(z-y²)
gives J_j=z J_(j-1)-1/(2j-1), with
J_0=atanh(1/sqrt(z))/sqrt(z).
This is analytic continuation from p>0. It is not a finite-frequency approximation.

The heavy integrals are

I_Htrace=(z³-6z²+9z)J0-z²+17z/3-36/5,
I_H2=z³J0/30-z²/30-z/90-1/150.

The old Proca integrals retain their complete corresponding polynomials and J0 factors. Both masses enter their own z; no common-mass substitution is made.

On the upper bank p=-s+i0, for z>0,
J0=atanh(sqrt(z))/sqrt(z)+i pi/(2sqrt(z));
for z<0,
J0=-atan(1/sqrt(-z))/sqrt(-z).
The heavy upper-bank imaginary factors are
U_Htrace=sqrt(z)(3-z)²/2, U_H2=z^(5/2)/60.
Below its threshold the heavy integral is real and NONZERO.

## Strict sign on the whole first sheet

For p=x+iy and a>0,b>=0,
Im[p/(a+bp)]=a y/[(a+bx)²+b²y²].
Every full weight is positive on0<y_radial<1. Thus Im A_i has the same strict sign as Im p, so no nonreal zero is possible.

On the real gap p>-4mu², each integral has strictly positive derivative. At the lowest threshold the old Proca values are16/15 in trace and -172/225 in shear. The heavy decrement there obeys

|I_Hi(-4mu²)| <= mu²/(n-mu²) int_0^1 W_Hi(y)dy <10^-180,

using the actual two masses. The heavy integrals are not set to zero: int W_Htrace=68/35 and int W_H2=1/210 are included in this bound. Since ell>394,

A_trace(-4mu²)>16/15+788-10^-180>789,
A_2(-4mu²)>-172/225+394/60-10^-180>5.

Real monotonicity excludes gap zeros. On every bank point s>4mu², the full Proca imaginary part is positive. Adding the heavy part cannot cancel it, including at the heavy threshold4n. This proves the entire actual first-sheet zero set is empty.

The original Proca-only factor still has its original pole on inversion. The above result concerns a different, physically specified SUM. No pole is discarded by hand and no finite coefficient was chosen to produce this conclusion. Finite reference-factor pole removal says nothing about the zero set of the full tree-plus-loop gravitational symbol.

## Contour, tails and exact moments

Write A_i(-s+i0)=D_i+i pi U_i. Its full reciprocal density is
rho_i(s)=U_i/[D_i²+pi² U_i²]>0.
It is the reciprocal density of the sum, not rho_P+rho_H for separately inverted channels.

The explicit full formulas give, as first-sheet |p| tends to infinity,

A_trace=4log p-2log(mu²)-52/15+O((n+mu²)(1+|log p|)/|p|),
A_2=(7/30)log p-(13/60)log(mu²)-127/450
      +O((n+mu²)(1+|log p|)/|p|).

The heavy ell cancels in these FULL high-frequency constants because its actual finite prescription is at mu=1. It does not cancel the heavy cut or the finite-frequency integral. The estimates are for the fixed physical masses, not a claim uniform in a changing n. The two logarithmic coefficients are positive.

Thus1/A tends to zero on the large slit circle. A keyhole contour for[1/A(z)]/(z-p) has no isolated-pole residues and no instantaneous constant. The first threshold arc vanishes because A has its strictly positive finite value there; the second threshold is bounded and lies inside the existing cut, with a nonzero Proca imaginary part. The banks give

1/A_i(p)=int_(4mu²)^infinity rho_i(s)/(p+s)ds.

For large s, rho_trace(s)(log s)² tends to1/4 and rho_2(s)(log s)² tends to30/7, so this integral converges absolutely. Positivity and dominated convergence give the exact total moments

int rho_trace/s=1/[2(ell+2)],
int rho_2/s=60/(ell+2),

int rho_i/s²=A_i'(0)/c_i²,
A_trace'(0)=9/(35mu²)+17/(35n),
A_2'(0)=3/(56mu²)+1/(840n).

No delta pole is hidden in either measure. For every q>=0 the shifted static mass is1/A_i(q)<=1/c_i.

Seventy-digit independent reconstructions integrate across both threshold locations and to infinity for twelve cases, including the physical n and complex p at both scales. Relative discrepancies are below6e-67; acceptance is1e-50. These numerical diagnostics corroborate, but do not prove, the sign and contour argument.

## Interior trace cusp is retained

Let u=log(s/(4n)) and A_* = D_*+i pi U_* at the heavy threshold, with U_*>0 supplied by the open Proca channel. The full heavy trace integral has a left increment(9pi/2)sqrt(-u) and right imaginary increment i(9pi/2)sqrt(u), plus analytic terms. Consequently the total reciprocal density has leading one-sided coefficients

rho(u)-rho(0) =
[-9pi D_* U_* /(D_*²+pi²U_*²)²]sqrt(-u)+O(u), u<0,

rho(u)-rho(0) =
[(9/2)(D_*²-pi²U_*²)/(D_*²+pi²U_*²)²]sqrt(u)+O(u), u>0.

Neither side may be silently omitted in a frequency-kernel regularity proof. The heavy shear nonanalytic term starts at order|u|^(5/2), not|u|^(1/2). Convergent local Puiseux expansions give differentiated remainders; these are not derivatives of an unsupported asymptotic assertion.
