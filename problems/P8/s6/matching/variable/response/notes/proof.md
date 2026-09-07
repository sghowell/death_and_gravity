# Proof of fixed-slice physical tensor transfer

All matrices use the Euclidean subordinate norm. We repeatedly use the sum
of absolute entries as an upper bound; no unsupported entrywise-to-operator
identification is made. All time derivatives are `d/du`, `u=T/tau`. The
model, physical metric, exact background and scalar source are those of
frozen VARIABLE; this proof changes neither the action nor its preparation.

## 1. Full operator and continuous coefficient enclosures

Write `d=1+u^2`, `c=2+delta`, and retain the parent's exact canonical fields
and inverse map. Define

\[
l''+A l+C Q+u d_c Q'=j_l,\qquad
Q''+B Q+E l+u f_c l'=j_Q,
\]
\[
A=q^2c_L^2-N_\Sigma,\quad
C=q^2D-2\omega'-2\omega\theta_\Sigma,\quad d_c=-2\omega/u,
\]
\[
E=q^2D-2\omega\theta_\Sigma,\quad f_c=2\omega/u,\quad
B=\mathfrak m^2+b_{\rm an},\quad
b_{\rm an}=q^2c_H^2-N_R-4\omega^2.
\tag{6}
\]

`d_c,f_c` have removable values at zero. `b_an` is analytic on the whole
closed box `|u|<=.1`, `2<=c<=2.01`, `1<=kbar<=2`. The pole-subtracted
quantity used below is a *different* function,

\[
b_\delta=B-\frac{80}{\delta+8u^2}.
\tag{7}
\]

It is bounded, but need not extend continuously to `(u,delta)=(0,0)`.
Indeed the mass-only remainder tends to `-32` on `u=0`, and `-141` on
`delta=0,u->0`. We never bound its delta derivative on the shrinking chart.

Here are uniform bounds, with their derivation rather than point samples:

\[
|A|,|C|<28,\quad |d_c|,|f_c|<14,\quad
|E|\le160(\delta+u^2),\quad |b_\delta|<200.
\tag{8}
\]

The E inequality is strict when `delta+u^2>0`; its analytic origin has
both sides zero. The pole-containing expressions are used only off that
origin, or in their specified bounded estimates.

The parent's exact margin gives `d^12<8/7`. Let
`R=(cd^12-8)/(cd^12+8)`. Then `|R|<=1`, and
`R'=384cu d^11/(cd^12+8)^2`, whose absolute value is at most `128/105`
using `c/(c+8)^2<=1/36` on `[2,4]`. Hence
`|theta_Sigma|<=6|u|`, `|theta_Sigma'|<7`, and `|N_Sigma|,|N_R|<8`.
The increasing function `sqrt(2c)/(c+8)` is bounded by `sqrt(8)/12`
on this interval. It follows that `|omega/u|<48/7<7`. Its logarithmic
derivative in absolute value is less than `20|u|`, so `|omega'|<9`.
Both squared gradient diagonals are convex combinations of `1` and
`c^2d^8/4<32/7<5`. Also `sqrt(w1w2)<=1/2`, `q^2<=4`.
These establish the first two bounds and `|b_an|<30`.

For the vanishing coupling, the two rectangle derivatives of `c^2d^8/4`
are bounded by `16/7` in delta and `256/7` in `v=u^2`, respectively. Thus
`c_f^2-1<37(delta+v)`. It follows that `|q^2D|<74(delta+v)` and
`2|omega theta|<=84v`, giving the `160` bound in (8).

For the pole set `S=delta+8v`, `D_m=(2+delta)(1+v)^4-2` and

\[
N=16(1-v)\left[\frac8{(2+\delta)(1+v)^{14}}+
                         \frac1{(1+v)^2}\right],\qquad
\mathfrak m^2=N/D_m.
\]

On `0<=v<=.01`, `N<=80`, `|N_delta|<=32`, `|N_v|<=1008`, so
`|N-80|/S<=126`. The exact positive polynomial difference satisfies

\[
0\le D_m-S\le(41/10)\delta v+(121/10)v^2.
\]

Use `delta v/S^2<=1/32` and `v^2/S^2<=1/64` to obtain

\[
|\mathfrak m^2-80/S|\le126+
80\left[\frac{41}{320}+\frac{121}{640}\right]
=1211/8<152.
\tag{9}
\]

Adding `|b_an|<30` proves (7)'s bound. Away from zero, `r=|u|<=.1`,
separate differentiation of this mass, not of (7), gives

\[
|\partial_\delta\mathfrak m^2|
\le\frac{32}{8r^2}+\frac{160}{64r^4}<\frac3{r^4}.
\tag{10}
\]

The remaining six analytic c derivatives of
`(A,C,d_c,f_c,E,b_an)` have absolute value `<8`. The verifier independently
reconstructs all of them by **bivariate interval Taylor arithmetic** from
the normalization formulas. Forty adjoining closed u intervals cover the
whole box, with full c and k intervals, not point evaluations. All interval
operations are exact Fractions; square-root endpoints follow integer
square roots, and final displayed bounds are rounded only outward. The
mixed `(u,c)` coefficient computes the c derivative of `omega'` and
`theta'`; no derivative is dropped. This independently verifies the source
projection and its c derivative are `<1` as well. Primary exact bridges
and the separately authored audit check the literal action/map formulas.

## 2. Two outer regular-singular systems

On either side `u=sigma r`, use the formulation's `Y` and `s=log r`.
Direct chain rule gives

\[
l_s=\sigma r l',\quad
(l')_s=-\sigma r A l
-\sigma r^{3/2}[(C+d_c/2)Y_H+d_c\mu P_H],
\]
\[
(Y_H)_s=\mu P_H,\quad
(P_H)_s=-\mu Y_H-
\frac{r^2(B-10/r^2)}\mu Y_H
-\frac{r^{3/2}E}\mu l-\frac{\sigma r^{5/2}f_c}\mu l'.
\tag{11}
\]

At delta zero the remainder in (11), beyond `A_0`, has entry sum at most

\[
r\{29+(35+14\mu)\sqrt r+(200/\mu)r
 +(160/\mu)r^{5/2}+(14/\mu)r^{3/2}\}<64r.
\tag{12}
\]

For `r<=.1`, use `3<mu<13/4`, `sqrt(.1)<1/3` to bound the braces by
`377/6<64`. This constructs the two full wave operators by the absolutely
convergent Volterra series (1). Both it and its inverse have norm at most
`exp(64r)`, and their difference from `r^A0`, respectively its inverse,
has norm at most `exp(64r)-1`. Truncating after degree n gives an error at
most `exp(64r)(64r)^(n+1)/(n+1)!`, an explicit computable interface.

We use the elementary bound, proved by `n!>=2*3^(n-2)` for `n>=2`,

\[
e^x\le1+x+\frac{x^2}{2(1-x/3)},\quad 0\le x\le2/3;
\qquad e^{2/3}\le41/21<2.
\tag{13}
\]

Consequently both outer wave norms are `<2` at `r<=a_*=.01`, and `<4`
at `r<=2a_*`. At `m=delta^(1/3)<=.001` their errors from the rotations
are `<=128m`.

For finite delta the difference of generators in (11) has norm bounded by

\[
\|M_\delta-M_0\|\le2\delta/r^2.
\tag{14}
\]

In fact (10) contributes `<delta/r^2`, while the remaining entries divided
by `delta/r^2` are bounded by
`8r^3+(12+8mu+8/mu)r^(7/2)+(8/mu)r^(9/2)+(8/mu)r^4`.
The total exact enlargement is `22993/22500<2`. No joint continuity of
`b_delta` is used. Conjugating by the delta-zero outer waves gives a
generator of norm at most `8delta/r^2`. Its integral over `[m,a_*]` in
`dlog r` is at most `4delta/m^2=4m`. Thus either orientation's finite-delta
outer coefficient transfer `L_+`, `L_-` satisfies

\[
\|L_\pm\|<2,\qquad\|L_\pm-I\|\le8m.
\tag{15}
\]

## 3. Exact inner transform, including both light components

Set `epsilon=sqrt(delta)`, `r_e=sqrt(u^2+delta/8)`,
`z=asinh(sqrt8*u/epsilon)` and `Q=sqrt(r_e) psi`.
Then `du/dz=r_e`. In `Z_i=(l,l',psi,psi_z/mu)` the universal generator
is zero on the light components and is exactly

\[
\psi_{zz}+[\mu^2+(3/4)\operatorname{sech}^2z]\psi=0
\tag{16}
\]

on the heavy pair. The exact full-minus-universal generator has entry sum
bounded by

\[
r_e\{29+(35+14\mu)\sqrt{r_e}+(200/\mu)r_e
+(1280/\mu)r_e^{5/2}+(14/\mu)r_e^{3/2}\}<72r_e.
\tag{17}
\]

For example the exact light-to-heavy forcing is `-r_e^(3/2) E l/mu`,
not a discarded term. Here `delta+u^2<=8r_e^2`; on `|u|<=m` the parameters
give `r_e<1/8`, `sqrt(r_e)<3/8`. The braces in (17) are bounded by
`6743/96<72`. Since `r_e dz=du`, their total integral is at most `144m`,
even though the length in z diverges logarithmically.

The propagator and inverse of (16) between *any* two real z values have
norm `<2`: remove its constant skew rotation and integrate the remaining
entry to get at most `3/(2mu)<1/2`. Conjugating the full perturbation
therefore has integrated norm at most `576m<=.576<2/3`. Variation of
constants and (13) first give a full propagator norm `<4` on every
subinterval. The **direct** Duhamel comparison then gives
`||F_i-F_PT||<=4*2*integral72r_e dz<=4*2*144m=1152m`:

\[
\|F_i\|<4,\qquad \|F_i-F_{PT}\|\le1152m.
\tag{18}
\]

This last bound is not inferred from the generally weaker exponential
expression `2(exp(576m)-1)`.

At `u=sigma m`, put `eta=delta/(8m^2)=m/8`. The exact overlap map is

\[
\psi=(1+\eta)^{-1/4}Y_H,\qquad
\psi_z/\mu=\sigma\{(1+\eta)^{1/4}P_H+
\eta Y_H/[2\mu(1+\eta)^{3/4}]\}.
\tag{19}
\]

It is the identity on `(l,l')`. Relative to
`J_sigma=diag(1,1,1,sigma)`, both this map and its inverse differ by
norm at most `eta=m/8`, by the slopes `1/4` of the fractional powers and
`mu>3`; both actual norms are `<2`. Thus both charts overlap at the
specified matching radius without replacing physical derivatives.

## 4. Exact connection, radial orientation and controlled tails

For (16), take the right Jost solution

\[
e^{i\mu z}\,{}_2F_1(-1/2,3/2;1-i\mu;(1-\tanh z)/2).
\]

The hypergeometric argument stays in `(0,1)`. The hypergeometric differential
equation directly verifies (16). The z=0/z=1 connection formula gives
the `A,B` in (3) at the left end. Although the two numerator parameters
differ by an integer, this is not the infinity pair to which that issue
applies; the local exponent difference used here is `-i mu`, nonintegral.
Gamma recurrence and reflection give
`B=i/sinh(pi mu)` and `|A|^2=coth(pi mu)^2`. The associated Ferrers
representation, if used, is the real-argument one, not the x>1 Legendre
chart. This is a classical scattering normalization, not a quantum state.

As `|x|=|u|/epsilon->infinity`, `z=sign(u)log(4sqrt2 |x|)+o(1)`.
After normalizing the right power to `x^(1/2+i mu)`, its left asymptotic
is `B |x|^(1/2+i mu)+A(4sqrt2)^(-2i mu)|x|^(1/2-i mu)`.
Changing from x powers to physical u powers and then inverting the
left-from-right matrix gives exactly (3). The same positive-frequency
phase has opposite radial exponents on the two sides. The unitary map
`(C1,C2)->((C1-i C2)/sqrt2,(C1+i C2)/sqrt2)` identifies these powers with
the outer real coefficient norm. Both phases in (3) are retained.

Here is an explicit tail bound. Let `Z=asinh(sqrt8*m/epsilon)`.
Then `1-tanh Z<=eta/2`, so each normalized Jost tail has integrated
perturbation at most `eta/8`; each wave matrix and inverse differs from
its limiting rotation by at most `eta/4`. Replacing both tails changes
the finite PT transfer in rotated coordinates by at most `2eta`.
Also

\[
0\le Z-\log(4\sqrt2 m/\epsilon)
=\log[(1+\sqrt{1+\eta})/2]\le\eta/4.
\]

Replacing the two phase rotations therefore costs at most `mu eta`.
Together these errors are at most `(2+13/4)eta<2m`.

The energy bound also gives `||S||<2`. Reflection is nonzero but small:
`pi mu>9` and the exact Taylor sum `sum_0^8 9^n/n!>2002` imply
`0<csch(pi mu)<1/1000`. Even with reflection set to zero, the nonzero
transmission amplitude retains its log-epsilon phase. In particular no
fixed transfer limit is asserted in place of (2)'s delta-dependent S.

## 5. Full matched estimate and actual physical norms

At the matching radius the exact central coefficient matrix is

\[
M=\Psi_+(m)^{-1}B_+^{-1}F_iB_-\Psi_-(m).
\]

Compare to the same product with the outer wave matrices replaced by
rotations, `B_sigma` by `J_sigma`, and `F_i` by `F_PT`. Replacing one
factor at a time, using ideal factors on its left and actual factors on
its right, gives the following coefficient multiplying m:

\[
128(2)(4)(2)(2)+(1/8)(4)(2)(2)+1152(2)(2)
+(1/8)(2)(2)+128(2)=17925/2<9000.
\]

Surround this middle product by the two outer transfers (15). The central
error becomes at most `36000m`; replacing those outer factors costs at
most `48m`, and the connection/tail replacement at most `2m`. Thus
`36050m<40000m`, proving (2).

Every endpoint retains `C_delta`, not a silently substituted `C_0`.
For explicit norm equivalence, on `|u|<=2a_*` the rational margins give
`1/2<f_Sigma<2`, `1/4<f_R<1`, `|theta|<=6|u|`, and `|w_i'|<7|u|`.
The physical map's entry sum is less than 32; direct inverse
differentiation of `l=f_Sigma(w1 g+w2 f)`, `Q=f_R(f-g)` gives an entry
sum less than 10. At `r=a_*`, `||W||<40`, `||W^-1||<14`, and the two
outer wave norms are `<2`. Hence `||P||<2560<3000`, `||P^-1||<280`.
For raw proper derivatives their additional norm factors are respectively
`max(1,tau^-1)` and `max(1,tau)`; canonical variables normalized without
dividing by M additionally carry that explicitly reversible M factor.

This proves a physical four-Cauchy-data result. It does not identify any
two-dimensional subspace as an EFT automatically. The small-error example
`delta<=10^-21` is a finite-parameter example of (2); the upper box bound
40 at `delta=10^-9` is deliberately recorded as a non-small control.

## 6. Actual g-source loading and paired prepared data

The unit-polarization TT stress and source normalization are as stated in
the formulation. Its energy divergence is H times the spatial trace and
its spatial divergence is proportional to `k_j e^{ij}`, both zero;
arbitrary smooth time dependence therefore preserves linear conservation.
Literal variation of `+a^3 sigma gamma_g/2` gives the two nonzero canonical
projections. In particular the heavy source is not removed by `omega(0)=0`.

On the fixed punctured interval `[-2a_*,-a_*]`, set delta=0 only in its
regular coefficients. For a desired outer coefficient C at `-a_*`, solve
the full homogeneous operator backwards there using (1). Choose a fixed
smooth cutoff zeta, zero near `-2a_*` and one near `-a_*`; for definiteness
take `zeta(u)=s(2(u+2a_*)/a_*-1/2)`, where `s(v)=0` for `v<=0`, `s(v)=1`
for `v>=1`, and `s(v)=exp(-1/v)/(exp(-1/v)+exp(-1/(1-v)))` otherwise.
Define `gamma_f=zeta gamma_f,target`, and use the literal positive spring
`U_TT=K_R mathfrak m^2` to set

\[
\gamma_g=\gamma_f+
[(K_f\gamma_f')'+G_f\bar k^2\gamma_f]/U_{TT},
\]
\[
\sigma_0=\frac4{a^3}
[(K_g\gamma_g')'+G_g\bar k^2\gamma_g+U_{TT}(\gamma_g-\gamma_f)].
\tag{20}
\]

Here `K_g=a^3/8`, `K_f=b^3/(8c)`, `G_g=a/8`, `G_f=cb/8`.
The unsourced f equation vanishes identically by construction. The g
equation equals its correctly normalized source. Near the left both
fields vanish; near the right they equal the homogeneous target because
that target satisfies the f equation as well. Thus `sigma_0` is a fixed
real smooth compact source with finite exact norm `S=integral|sigma_0|du`.
Its definition uses no delta-dependent reset or pulse compression. Volterra
iteration with the factorial tail, the explicit smooth cutoff and regular
coefficient derivatives provide a computable source functional; no numeric
value of S is invented or omitted from the bound.

Use that same source for every admitted positive delta. In outer normalized
coordinates on `[a_*,2a_*]`, both wave norms are `<4`. The coefficient
perturbation is bounded by `32delta/r^2`, with integrated norm at most
`160000delta`. The source and its c derivative have Euclidean log-time
norm at most `2r|sigma_0|`; after conjugation the bounds are `8r|sigma_0|`
and `8delta r|sigma_0|`. Since `|dlog r|=|du|/r`, the exact finite-delta
retarded coefficient is bounded by `16S`. Its difference from delta zero
is at most

\[
[8+160000(16)]\delta S=2560008\delta S<3000000\delta S.
\]

Combining this with (2) and `||S_epsilon||<2` gives (5). The specified
targets `C=(0,0,1,0)` and `C=(1,0,0,0)` therefore give a paired incoming
relative-mode and prepared regular-light control. The latter has a
vanishing bounded leakage, but the estimate does **not** prove a sharp
`O(delta)` leakage rate or an omission error below the locked-cone excess.

For a nonzero relative target, the off-diagonal transmission phase in
(3) yields nonconvergent output along two sequences of delta tending to
zero with phases separated by pi. Their ideal heavy output separation is
`2|A| ||C_H||`; (2) and (4)'s errors vanish. This is a full physical data
statement because `P_delta^+` tends to an invertible punctured endpoint
map. It is also observable using g alone after the source switches off:
the g equation reads `g''=alpha g+beta g'+j f`, `j=U_TT/K_g>0`.
Differentiating once gives a map from `(g,g',f,f')` to `(g,g',g'',g''')`
with determinant `j^2>0`, finite and nonzero at the fixed endpoint as
delta tends to zero. No assertion about only the pair `(g,g')` or a
low temporal-band measurement is substituted for this exact finite-jet map.

## 7. What the result does not decide

The exact source response is infinitesimal: source-amplitude differentiation
precedes the delta limit. No tensor backreaction or quantum production is
calculated. A fixed-duration compact source is not a low temporal-frequency
band. The fixed spatial momentum band does not establish a temporal gap,
and `delta=0` is not an admissible full action through the center. No
uniform complete-background or light-only source/state matching estimate
is inferred from the full high cone or from the nonconvergent heavy phase.
In particular a prepared-light family can have suppressed relative response;
the sharp suppression order and its low-energy significance remain open.
This bounded response theorem is not an original-S6 UV or positivity verdict.
