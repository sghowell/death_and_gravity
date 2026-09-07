# Complete coefficient-space and physical-response proof

All primes in this note mean `d/du`, not proper-time differentiation.
The field and source normalizations are those of frozen S6.21. In
particular canonical fields are divided by `M`, and a physical derivative
is `tau^-1` times its displayed u derivative. Constants below are exact
rationals. Decimal approximations, where given for orientation, are not
inputs to any inequality.

## 1. Literal coupled operator and analytic pole removal

Put `d=1+u^2,c=2+delta,K=(tau*kcom)^2`, and

```
D=c*d^4-2,
N=16(1-u^2)[8/(c*d^14)+1/d^2],
m^2=N/D.
```

The exact canonical equations from the physical action are

```
l''+A*l+C*Q+u*d_cross*Q'=0,
Q''+(m^2+b_an)*Q+E*l+u*f_cross*l'=0.                 (1)
```

Every coefficient other than `m^2` is jointly analytic near the origin.
In the notation of the parent,

```
theta=6u(c*d^12-8)/(d(c*d^12+8)),
omega=-24sqrt(2c)*u*d^5/(c*d^12+8),
w1=c*d^12/(c*d^12+8), w2=8/(c*d^12+8),
c_f^2=c^2*d^8/4, sqrtw=sqrt(8c)*d^6/(c*d^12+8),
q_spatial^2=K/d^4,
A=q_spatial^2(w1+c_f^2*w2)-theta'-theta^2,
C=q_spatial^2*sqrtw*(c_f^2-1)-2omega'-2omega*theta,
E=q_spatial^2*sqrtw*(c_f^2-1)-2omega*theta,
d_cross=-2omega/u, f_cross=2omega/u,
b_an=q_spatial^2(w2+c_f^2*w1)+theta'-theta^2-4omega^2.
```

`b_an` is not the merely bounded, path-dependent remainder obtained by
subtracting `80/(delta+8u^2)` from the full heavy coefficient.

For actual positive delta, `D>0`. Substitute the exact identity `Q=Dq`:

```
l''+A*l+(C*D+u*d_cross*D')q+u*d_cross*D*q'=0,
Dq''+2D'q'+(D''+N+D*b_an)q+E*l+u*f_cross*l'=0.       (2)
```

These equations, unlike the original mass coefficient, are analytic on
the complex proof polydisc. Solving (2) there and restricting to actual
positive delta gives exact solutions of (1). There is no inference that
the delta=0 full action exists.

## 2. Analytic inverse, derivative compositions and full contraction

Let `R=1/20,v=R^2=1/400`. For a series
`f(delta,u;K)=sum f[j,n](K)delta^j u^n`, use the Banach algebra norm

```
||f|| = sum_(j,n>=0) |f[j,n](K)| v^j R^n,          (3)
```

uniformly for complex `|K|<=4`. Equivalently the unknown is a holomorphic
K-valued map into this coefficient-l1 Banach space, bounded in its norm
uniformly on that disk. Multiplication is submultiplicative. Every norm
below is a coefficient norm, not merely a pointwise real-axis supremum.

The leading operator is

```
L0=(delta+8u^2)partial_u^2+32u partial_u+96.
```

Its diagonal on `delta^j u^n` is `d_n=8(n^2+3n+12)`. After diagonal
inversion the remaining delta shift has weight ratio
`(n+2)(n+1)/[8(n^2+3n+12)]<1/8`. The equality `v=R^2` is essential.
Thus its inverse exists on (3), with norm at most `(8/7)/96=1/84`.
On a polynomial monomial the shift chain is finite, which gives the exact
triangular replay in `analytic.inverse_monomial`. On the complete space,
the Neumann series proves convergence and uniqueness.

Needed compositions, proved term by term before the same Neumann inverse,
are

```
||L0^-1(Pq)|| <= ||P|| ||q||/84,
||L0^-1(Pq')|| <= ||P|| ||q||/(70R), if P has u-degree>=1,
||L0^-1(Pq'')|| <= ||P|| ||q||/(7R^2), if P has u-degree>=2.
```

For the second line use `n/d_(n+m-1)<=1/80` for integer `n>=0,m>=1`;
the maximum at `m=1,n=3,4` is `1/80`. For the third use
`n(n-1)/d_(n+m-2)<=1/8` for `m>=2`. These estimates are not a claim that
uncomposed differentiation is bounded on the original radius.
Twice integrating from zero, denoted `J^2`, satisfies

```
||J^2 f||<=R^2||f||/2,
||J^2(Pq')||<=R||P|| ||q||/6  if P has u-degree>=1,
```

because `n/[(n+m)(n+m+1)]<=1/6` for the relevant integers. All these
inequalities include their zero-index endpoints.

For transparency the continuous coefficient majorants are specified
completely here and replayed using both SymPy rationals and independent
Fractions. Write `cp=2+v,cm=2-v,dp=1+v,dm=1-v`,

```
F=cp*dp^12-2, Z=10-F, Fv=12cp*dp^11,
T=6(6+F)/(dm*Z),
Tv=6[Fv/(dm*Z)+(6+F)/(dm^2*Z)+(6+F)Fv/(dm*Z^2)],
O=48dp^5/[(1-v/2)Z],
Ov=48/(1-v/2)[5dp^4/Z+dp^5 Fv/Z^2],
Ntheta=T+2v*Tv+v*T^2, Oprime=O+2v*Ov,
W1=cp*dp^12/Z, W2=8/Z, Cf=cp^2*dp^8/4,
Wsqrt=4dp^6/[(1-v/2)Z], Qs=4/dm^4.
```

These follow by geometric series about the nonzero constant denominators.
For example `||(1+u^2)^-p||=(1-v)^-p`, and the denominator
`c*d^12+8=10+(c*d^12-2)` has reciprocal norm at most `1/Z`.
The binomial square-root coefficients are majorized by a geometric
series. `c_f^2-1` has zero constant term and positive remaining
coefficients, so its norm is `Cf-1`, not `Cf+1`. It follows that

```
||A|| <= Qs(W1+Cf W2)+Ntheta <9,
||C|| <= Qs Wsqrt(Cf-1)+2Oprime+2vOT <12,
||E|| <= Qs Wsqrt(Cf-1)+2vOT <60v,
||b_an|| <= Qs(W2+Cf W1)+Ntheta+4vO^2 <9,
||d_cross||,||f_cross|| <=2O<10.
```

Set

```
Dnorm=cp*dp^4-2,
UDprime=8v*cp*dp^3,
DppError=8cp*dp^3+48v*cp*dp^2-16,
Nerror=16[(1+v)(8/(cm*dm^14)+1/dm^2)-5],
tq=(Dnorm-9v)/(7v)+(16cp*dp^3-32)/70
    +(DppError+Nerror+9Dnorm)/84,
a=9v/2, b=(23/3)Dnorm+5UDprime, c=(1+5v)/7.        (4)
```

Here `||D-delta-8u^2||=Dnorm-9v`, and that difference has u-degree at
least two. Also `2D'-32u` has u-degree at least one. The alternating
Taylor coefficients of `N` give the displayed `||N-80||<=Nerror`.

The full fixed point is

```
l=l_init-J^2[A*l+(C*D+u*d_cross*D')q+u*d_cross*D*q'],
q=-L0^-1[(D-delta-8u^2)q''+(2D'-32u)q'
         +(D''+N+D*b_an-96)q+E*l+u*f_cross*l'].
```

In the pair norm `||l||+v||q||`, its two column sums are
`a+vc<1/2` and `b+tq<1/2` (approximately `.01162,.42331`). Therefore it
is a strict contraction on the entire Banach space. The fixed point
exists for every pair of light initial constants, depends holomorphically
on K, and solves (2) coefficientwise and classically in the open polydisc.
There is no free jointly analytic homogeneous heavy mode. Uniqueness is
within this prepared analytic class, not among all actual Cauchy data.

The coefficients preserve u parity. For the even solution `l_init=1`,
separating the initial constant and using the contraction gives

```
||l_e-1|| <=2(a+5v^2/7)<10v,
||q_e|| <=[5v/7+c*2(a+5v^2/7)]/(1-tq)<3v.
```

For the odd solution `l_init=u`,

```
||l_o-u|| <=2(a+vc)R<R/40,
||q_o|| <=c[R+2(a+vc)R]/(1-tq)<R/5.                (5)
```

The initial conditions hold for every delta; hence `l_e-1` has u-degree
at least two and `l_o-u` has u-degree at least three. The weighted-degree
zero equation is `96q_e[0,0]=0`. Thus the even q has no joint constant
term; the odd q has u-degree at least one. In particular
`|q_e(delta,0)|<=3delta`, a stronger fact than its disk bound `3v`.

## 3. Leading coefficients and a uniform finite-delta center remainder

Direct expansion of the full coefficients, including the moving weights,
gives

```
A(0,0)=K+18/5, f_cross(0,0)=-48/5,
E=(2/5)K*delta+[(16/5)K-864/25]u^2+weighted-degree>=4.
```

Solving the exact leading triangular system yields

```
l_e=1-(K+18/5)u^2/2+...,
q_o=(3/40)u+...,
q_e=-K[(4/55)u^2+(7/2640)delta]+... .             (6)
```

Thus under `delta=epsilon^2,u=epsilon*x`, the leading relative fields are

```
Q_o=epsilon^3[(3/5)x^3+(3/40)x]+...,
Q_e=-epsilon^4*K[(32/55)x^4+(31/330)x^2+7/2640]+... .
```

They solve the exact leading inner forced operator
`partial_x^2+80/(1+8x^2)`, but convergence comes from section 2, not from
these polynomial identities alone.

There is also a direct finite-domain bound for these full-field
remainders. Put `t=max(sqrt(delta/v),|u|/R)<=1`. Odd q after its displayed
term has weighted degree at least three; even q after (6) has degree at
least four. Since `D-delta-8u^2` has degree at least four, the same
coefficient-l1 norm proves

```
|Q_o-(delta+8u^2)(3u/40)| <=Dnorm*(R/5)*t^5,
|Q_e+(delta+8u^2)K[(4/55)u^2+(7/2640)delta]|
 <=3v*Dnorm*t^6.
```

These are uniform bounds on the stated actual slab; no unspecified
O-symbol is needed to use the inner polynomials.

At `K=0`, constant equal physical fields are an exact full solution:
`Q=0,l=f_sum(u,delta)/f_sum(0,delta)`. Its analyticity and uniqueness
identify it with the even branch. Consequently `q_e(K=0)=0` exactly.
Apply scalar Schwarz to every norm-one bounded linear functional on the
coefficient-l1 Banach space (3), then its dual norm formula. Since
`K -> q_e(K)` is Banach-holomorphic, vanishes at zero and has norm at most
`3v` on `|K|<=4`,

```
||q_e/K|| <=3v/4.
```

This is a bound in the **coefficient-l1 space**, not just a pointwise
supremum in delta. At u=0 the coefficients of delta^0 and delta^1 are
respectively zero and `-7/2640`. Therefore the sum of all delta^j terms
with `j>=2` obeys

```
|q_e(delta,0)/K+(7/2640)delta|
 <= (delta/v)^2 ||q_e/K|| <=300delta^2.             (7)
```

No `1/(1-delta/v)` factor is needed: the total absolute coefficient sum
is already bounded. A plain supremum-norm Taylor estimate would not give
this exact constant by the same argument. Equation (7) includes K=0 by
holomorphic continuation. For actual `K>0`, it also shows `q_e<0`.

The actual physical g equation at the center is
`g''+K g-mu_TT Q/f_R=0`, with
`mu_TT=128/[delta(delta+2)]`. On the even branch its exact multiplier is

```
G(0)/K=1-B*q_e(delta,0)/K,
B=32(delta+10)/[(delta+2)(sqrt(1+delta/2)-2delta*q_e(delta,0))]. (8)
```

For `0<delta<=10^-9,0<=K<=4`, `0<B<=160` and `|B-160|<=105delta`.
For example (7) gives `-q_e<=4[(7/2640)delta+300delta^2]`. Together with
`1<=sqrt(1+delta/2)<=1+delta/4`, a direct positive-denominator calculation
bounds `(160-B)/delta` by

```
104+(20+1280a0)delta+(640a0+384000)delta^2+192000delta^3 <105,
a0=7/2640.
```

All omitted denominator factors have been enlarged in the safe direction.
Combining (7),(8),

```
|G(0)/K-1-(14/33)delta|
 <=[48000+105*(7/2640)]delta^2 <48001delta^2.       (9)
```

This is a uniform finite-parameter statement for the selected branch, not
a merely formal coefficient. At `delta=10^-9` the displayed remainder is
less than `48001*10^-18`. The locked coefficient is exactly
`1+(4/5)delta+(6/5)delta^2/(10+delta)`, so its different first-order
coefficient is a genuine source-free preparation comparison. It does not
identify a high-frequency characteristic or a local arbitrary-source EFT.

## 4. Conserved physical symplectic normalization and the closed g equation

After the parent's boundary integration the canonical momenta are
`pi_l=l'-2omega Q,pi_Q=Q'`. The exact full conserved bilinear form is

```
Omega=l_e(l_o'-2omega Q_o)-(l_e'-2omega Q_e)l_o
      +Q_e Q_o'-Q_e'Q_o.
```

At zero, parity and the fixed light data imply
`Omega=1+delta^2 q_e(delta,0) q_o'(delta,0)`.
Section 2 gives `|q_e(delta,0)|<=3delta` and `|q_o'(delta,0)|<=1/5`.
Hence `|Omega-1|<=3delta^3/5`, in particular `9/10<Omega<11/10`.
The factor delta^3 uses the proved vanishing of q_e at the joint origin;
it cannot follow from its disk norm alone.

Write the actual physical map as `g=a_g*l+b_g*Q`,
`a_g=1/f_sum,b_g=-w2/f_R`. On `|u|<=1/100`, elementary positive
denominator estimates of the pinned map give

```
4/5<a_g<1, |b_g|<4, |a_g'|<=6/100, |b_g'|<=52/100,
D<=1/1000, |D'|<=1/5.
```

With `rho=(1/100)/R=1/5`, the coefficient bounds (5) and their stated
minimum u degrees give

| quantity | even upper bound | odd upper bound |
| --- | ---: | ---: |
| light deviation from `1,u` | `1/1000` | `1/100000` |
| derivative deviation from `0,1` | `1/5` | `3/1000` |
| q | `3/10000` | `1/500` |
| q' | `3/50` | `1/5` |
| Q | `3/10^7` | `1/500000` |
| Q' | `3/25000` | `3/5000` |

For example the derivative bound uses
`max_(n>=2)n rho^(n-1)=2rho`; the delta-only part of q_e contributes
at most `delta/v<=rho^2`. Including every derivative in the physical map
then gives

```
g_e>79/100, g_o'>79/100,
|g_e'|<27/100, |g_o|<11/1000,
g_e<501/500, g_o'<101/100.
```

Consequently the physical Wronskian
`W_g=g_e*g_o'-g_e'*g_o` has `3/5<W_g<2` on the full slab. Both physical
Cauchy components are retained, not just the field value at the center.
The exact source-free equation on this subspace is

```
g''+F*g'+G*g=0,
F=-W_g'/W_g,
G=(g_e'*g_o''-g_e''*g_o')/W_g.
```

Its kinetic normalization `K_eff=Omega/W_g` satisfies
`9/20<K_eff<11/6`, and `K_eff'/K_eff=F`. The action
`M^2/(2tau) integral du K_eff[(g')^2-Gg^2]` represents this fixed-K
source-free equation and its restricted symplectic form. At K=0 the
constant common mode implies `G=0` throughout the slab. Analyticity in K
therefore factors G by K locally along the real K interval. The
coefficient and kinetic function may depend on K: no spatially local
off-shell action, healthy unrestricted source response, or cone follows.

## 5. Endpoint inclusion and improved fixed-light error

Only this section imports the uniform transfer theorem from S6.21, whose
spatial band is retained as `1<=K<=4`. Let `r=1/100,mu=sqrt(39)/2` and
use its actual canonical-to-physical matrix `C_delta`, outer map `W`, and
full delta=0 wave matrix `Psi_side`. Define the prepared endpoint columns

```
H_delta^side=Psi_side(r)^-1 W_side(r)^-1 Z_delta(side*r),
Z_delta=(l,l',Q,Q') for the even/odd columns.
```

The corresponding physical endpoint matrix is
`P_delta^side=C_delta(side*r)W_side(r)Psi_side(r)`. Thus the actual
physical map cancels exactly in the definition of H. It is never replaced
by its delta=0 value. At delta=0 the analytic solutions have no relative
indicial data and the exact origin light normalization `(1,0),(0,1)`.
Uniqueness of the full outer Volterra solution gives
`H_0^side=L`, where `L` is the 4-by-2 inclusion `[I_2;0]`.

Coefficient-l1 parameter differences acquire a factor `delta/v`. At
either endpoint they imply the following bounds, each divided by delta:

| component | even | odd |
| --- | ---: | ---: |
| Delta l | `2/5` | `1/250` |
| Delta l' | `80` | `6/5` |
| Delta q | `3` | `4/5` |
| Delta q' | `24` | `80` |
| Delta Q | `9/2500` | `3/625` |
| Delta Q' | `74403/100000` | `3201/5000` |

The last two rows also use `|Delta D|<=2delta,|Delta D'|<=delta/10`
and the q and q' bounds at delta=0. The actual normalized outer vector is
`Y=(l,l',Q/sqrt(r),(uQ'-Q/2)/(mu sqrt(r)))`. Its component l1 norm
dominates the Euclidean norm, and the pinned inverse wave norm is below
two. The resulting even and odd endpoint column bounds are respectively
below `162delta` and `3delta`; their sum is below `200delta`. Thus

```
||H_delta^side-L||<=200delta.                      (10)
```

If `T_delta` denotes the full transfer in the actual P_delta endpoint
coordinates, exact homogeneous propagation gives
`T_delta H_delta^-=H_delta^+`. The frozen matched result gives
`||T_delta||<=2+40000delta^(1/3)<=42`. Equation (10) therefore proves

```
||T_delta L-L||<=200delta(1+42)=8600delta.          (11)
```

This is an O(delta) **upper rate** for fixed delta=0 regular-light data;
it supplies no sharp or nonzero leading leakage coefficient. The same
fixed old g-only source that targets `Lv` at delta=0 has an incoming
coefficient error at most `3*10^6 delta S`, where `S=integral|sigma_0|du`
is its actual norm. Hence its output differs from `Lv` by at most

```
delta[8600||v||+126000000 S].                      (12)
```

A small fixed-slice error is not a bound on every center derivative.
In particular (9) is not transferred to a fixed old source by (12).
No lower leakage rate, temporal-band estimate, or arbitrary-source EFT
claim is needed or inferred.

## 6. An actual finite, fixed-duration, delta-dependent g-source preparation

Let `(G,F)` be either exact analytic physical target, or a combination
with amplitude norm `V=|alpha_e|+|alpha_o|`. Define the fixed smooth cutoff

```
s(x)=f(x)/(f(x)+f(1-x)), f(x)=exp(-1/x) for x>0, f(x)=0 for x<=0,
zeta(u)=s(200(u+1/50)-1/2).
```

It is flat zero through `u=-7/400`, flat one from `u=-1/80`, with
transition strictly inside the fixed interval `(-1/50,-1/100)`.
Use `K_g=d^6/8,K_f=1/(c*d^6),G_g=d^2/8,G_f=c/(4d^2)` and
`U_TT=K_R*N/D`, where `K_R=d^6/(c*d^12+8)`. Set

```
f_loaded=zeta F,
g_loaded=f_loaded+[(K_f f_loaded')'+G_f K f_loaded]/U_TT,
sigma=4/a^3[(K_g g_loaded')'+G_g K g_loaded
            +U_TT(g_loaded-f_loaded)].           (13)
```

This solves the unsourced f equation and the physically normalized
g-source equation exactly, starts from zero, and equals the target after
the cutoff. The source is conserved because it has zero time components,
zero spatial trace and spatially transverse polarization, as independently
verified in the parent. The dictionary remains
`sigma=tau^2 Pi/M^2` with probe action
`(1/2) integral dT a^3 Pi gamma_g`.

Formula (13) must not be estimated by an apparent `delta^-1` spring.
Instead let

```
P=1/U_TT=D/(K_R N),
M_f=2K_f F'+K_f'F, N_f=K_f F,
A_s=P M_f, B_s=P N_f,
g_loaded=zeta G+zeta' A_s+zeta'' B_s.
```

`N` has constant term80 and `||N-80||<3`; thus P is analytic on the
entire proof polydisc, with `||P||<1/100`. Using both homogeneous target
equations cancels the spring exactly and gives
`sigma=(4/a^3)sum_(j=1)^4 S_j zeta^(j)`, where

```
S1=2K_g G'+K_g'G+(K_g A_s')'+G_g K A_s+M_f,
S2=K_g G+K_g A_s'+[K_g(A_s+B_s')]'+G_g K B_s+N_f,
S3=K_g(A_s+B_s')+(K_g B_s)',
S4=K_g B_s.                                      (14)
```

This identity is replayed with independent symbolic generic coefficients;
no target samples or division by its field value occur.

For explicit bounds shrink only the u coefficient radius to `R/2=1/40`,
keeping the delta radius v unchanged. A coefficientwise binomial estimate
gives

```
||partial_u^j f||_(R/2,v) <= j! 40^j ||f||_(R,v).
```

Indeed `binomial(n,j)t^(n-j)(1-t)^j<=1`, here `t=1/2`.
The rational coefficient bounds already give
`||K_g||,||K_f||,||G_g||<1`, `||P||<1/100`, and
`||G||,||F||<5V`. For the last assertion use
`||a_g||<2,||b_g||<4,||b_f||<1`, (5), and `||D||=Dnorm`.
Leibniz with the derivative bounds through order three bounds every S_j
in (14); `||4/a^3||<5`. The explicit integer/rational results are
`source.coefficient_bounds()`, independently reconstructed without SymPy.

To avoid an assumed smooth-cutoff norm, write
`f^(j)(x)=exp(-y)P_j(y)`, `y=1/x`, where
`P_0=1,P_(j+1)=y^2(P_j-P_j')`. The inequality
`y^n exp(-y)<=n!` gives, for `j=0..4`, the valid bump bounds
`(1,2,36,1584,129600)`. At least one of x and 1-x is at least one half;
on the transition the denominator is at least `exp(-2)>1/9`. The last
strict inequality follows from `e<3` by its factorial series. Put

```
I0=9,
In=9 sum_(j=1)^n binomial(n,j)*2*B_j*I_(n-j),
Z_n=200^n sum_(j=0)^n binomial(n,j)*B_j*I_(n-j), n=1..4.
```

Then `||zeta^(n)||_infinity<=Z_n`. Flat endpoints ensure all integration
and matching data are smooth. Denote the completely specified rational
source constant by `C_sigma=sum C_j Z_j`, with C_j the coefficient bounds
including the factor five. For every actual member,

```
||sigma_delta||_infinity<=C_sigma V,
integral|sigma_delta|du <=C_sigma V/100,
||sigma_delta-sigma_0||_infinity<=400delta C_sigma V. (15)
```

The last inequality again uses the coefficient-l1 norm with the unchanged
delta radius v: every nonconstant delta coefficient contributes at most
delta/v. The cutoff is not analytic in u, but it is fixed, and (14)
contains only its four bounded derivatives multiplying analytic
coefficients. This is why analytic coefficient estimates apply to (15).
The norm is finite and intentionally very coarse, not a claim of a
small finite-amplitude stress. The source is a linear probe evaluated at
zero amplitude before delta limits. Retuning it with delta is explicit;
a single fixed source is instead covered by (12).

## 7. Limits of the result

This proves a full-parent prepared analytic two-mode sector, not only an
inner toy polynomial. Its finite remainder, positive physical symplectic
normalization and actual source loading are independently testable.
It does not prove that arbitrary incoming data lie in that sector, that a
fixed source produces its special center coefficient, or that a quantum
vacuum chooses it. Compact temporal support does not imply a low temporal
band. No scalar/vector, nonlinear source, interacting-field, global
continuation, original C/D matching, UV-positivity or universal light-EFT
verdict follows. The frozen parent and all its source bytes remain
unchanged; the original S6 and P8 completion gates remain separate.
