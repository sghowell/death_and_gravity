# Whole original source and mixed real-time/complex-lapse bounds

## Exact binding, not a bounce jet

Write R=R_full(u,N), F=F_full(u,N), and j the unchanged mass-normalized
heavy source. U=R^(-3/4). The original primitive satisfies
I_N=3 R^(-7/4) R_u R_N/(4N), I(u,1)=0. The actual time coefficients are

B=-U R_u/(2N)-I,
Fhat=U[F+9 R_u^2/(16 R N^2)]-I_u/N,
lower=-3 Hclock(R-1)/N.

The whole normalized reduced Hamiltonian is N times

-3[p-B-(R-1)G]^2/(4 R^(1/4)) - Fhat
+ G^2/(2U) - lower G +2 sh R^(-1/4)
+[(1/10+dp)^2+ph^2]/(2U)
+R^(-1/4)(gm+gh)/2
+U[mu eta^2/2-j eta/10^100]
-R^(3/4) curv/2
+el R^(1/4)/2+ma R^(1/4)/4+R^(-1/4)wm/2.

The twelve variables are the actual complete density invariants, not
independent canonical oscillators. The definitions of sh,el,ma,wm,
gm,gh,curv retain the spatial metric, full lifted momentum, Gauss,
vector, both matter and all spatial curvature terms. eta=10^100 H,
mu=n/10^200 and the unperturbed matter density is1/10.

source.py literally substitutes these coefficients and all twelve
invariants into the unchanged S257 whole Hamiltonian and temporal
solution. The residuals vanish without replacing the original R,F
by a bare clock. The lower and I_u terms remain for nonzero real time.

## Uniform original coefficients

Use the S256 complete joint analytic source neighborhood. Set
r=1/10000, distance=11/10000. At real |u|<=T and complex
|N-1|<=rho, both center distances plus r are below distance.
The whole map N->N^-2 has image distance bounded by
distance*(2+distance)/(1-distance)^2<1/400.

The full original tree, fixed vacuum constants and all switch/heavy
terms give analytic moduli R<2, F_analytic<10^7,
j<10^-100 on that neighborhood. No term is removed merely because its
localizer is small. In particular the original constant and complement
prefactors10^301 and10^207 multiply the entire switch bound. For the
entire heavy localizer use exp(-x)<=8!/x^8 with the original coefficient
prefactor, giving the displayed bound below10^-2700.

For mixed derivatives i times real u and k times complex N, i+k<=5,
Cauchy multiplies the analytic bounds by i! k! r^(-i-k).
The original smooth fixed profiles are NOT asserted holomorphic in u.
Use their given real C5 derivative bound, with N-only Cauchy and the
full profile prefactor8, to add
8 PROFILE_BOUND k! r^-k to the F bound. All21 actual rational rows
give |partial_u^i partial_N^k R|, |...F|<10^30 and |...j|<1.
The proof needs no sixth profile derivative or complex-time extension.

The full clock identities are R(u,1)=1 and
R_N(u,1)=-2/(1+u^2)^3; the normalized source vanishes at N=1.
The entire source remainder and its first X jet vanish at X=1, checked
from the actual inherited source. Thus |R-1|<=10^30 rho<1/100 on
the lapse disk and the nonzero principal powers are well defined.
The actual rational Hclock=4u/(1+u^2) and its first3 derivatives
satisfy the separate interval bounds used below.

## Primitive and derivative majorant

Only derivatives of I containing at least one N derivative are replaced
by derivatives of the exact I_N formula. At fixed real u each pure
time derivative of I is the actual integral from1 to N of the
corresponding pure time derivative of I_N. On the straight complex
segment, its modulus is at most rho times the complete integrand
majorant. Setting I_u to zero off the clock would be incorrect.

branch.py differentiates the whole constraint C=partial_N Hbar and
bounds Hbar,C,C_N,C_u,C_NN,C_Nu,T and every C_z,C_Nz,C_uz,C_zz:
187 expressions, including every ordered second invariant contact.
The recursive engine uses triangle/product rules, rational coefficients,
|N|,|R|<2 and |N^-1|,|R^-1|<2. Fractional positive or negative
powers are rounded OUTWARD by the integer ceiling of their exponent.
Unproved derivative orders, reciprocals, fractional bases, inexact
coefficients and unknown atoms are rejected.

All187 rational majorants are<10^200, with the complete Hamiltonian
modulus<10^62. These bounds use the actual source jets, not an
unproved generic smoothness assumption. At u=0 the literal complete
constraint restricts to the unchanged S266 constraint, including its
R_uu and primitive-time contact; the symbolic identity is checked.

