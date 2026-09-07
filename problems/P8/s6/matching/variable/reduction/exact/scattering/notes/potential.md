# Uniform canonical own-f pole bound

This note proves a finite coefficient estimate for the exact own-f branch
of frozen S6.33, report SHA-256
`727811a8563584cc511e749ff038eadaa416f02aae268a15fa4f90eb8ca4b6b3`.
The physical metric and clock are the same external off-shell probe
`g=eta`, `theta=T`; this is not the full parent/CD background. The positive
branch is unchanged. Its `delta=0` extension is used only to prove uniform
bounds and limits, not as a literal action.

## 1. Canonical operator and the stated domain

Let `u=T/tau`, `v=u^2`, `c=2+delta`. On

\[
|u|\le1/100,\qquad0<\delta\le1/100,                                 \tag{1}
\]

the exact branch has `b>0`, `N>0`. Write `k=b^3/N`, and let `Q_tensor` denote
the hidden tensor amplitude; it is unrelated to the rational profile
numerator `Q_profile`. The physical metric tensor `q` is prescribed.
The own-f equation from the actual quadratic tensor action is

\[
(k Q_{{\rm tensor},u})_u+2b_1b(Q_{\rm tensor}-q)=0.
\]

With `Y=sqrt(k)*Q_tensor` (up to an irrelevant positive constant canonical
normalization), the exact equation is

\[
Y_{uu}+V_\delta Y=\sqrt{k}\,\frac A D q,\qquad
V_\delta=\frac A D-\frac{(\sqrt{k})_{uu}}{\sqrt{k}},\qquad
A=\frac{2Q_{\rm profile}N}{b^2},\quad b_1=Q_{\rm profile}/D.          \tag{2}
\]

Thus the first derivative is removed, but neither the canonical source
weight nor an independent homogeneous solution is discarded. The theorem
proved here is solely

\[
\boxed{\quad
\left|V_\delta(u)-\frac{16}{\delta+8u^2}\right|<44
\quad\text{throughout (1).}\quad}                                   \tag{3}
\]

`V_delta` is dimensionless in the `u` equation. In physical proper time its
potential is `V_delta(T/tau)/tau^2`: the pole is
`16/(tau^2*delta+8T^2)` and the remainder bound is `44/tau^2`. The constant
Einstein coefficient `M^2` cancels from this operator. This statement is
not a Green prescription, a source-transfer theorem, or a cutoff estimate.

## 2. Regularized background functions

Use the S6.33 definitions

\[
\begin{aligned}
d&=1+v,&J&=cd^4(1-7v)+12v,&Q_{\rm profile}&=32(1-v)/(cd^{14}),\\
R&=d^8J/[8c(1-v)],&D&=c-2/d^4,&L&=R_v/R,\\
H&=D(J_v/J-6/d)-D_v,&F&=\tfrac23(v\zeta H-L),\\
b&=[(1-vD\zeta)/R]^{1/3},&\mathcal T&=3bF^2/(2Q_{\rm profile}).
\end{aligned}                                                        \tag{4}
\]

The actual branch is `zeta=T(v,zeta,c)`. S6.33 proves its unique analytic
solution throughout the closed enclosure

\[
v\in[0,1/10000],\quad\delta\in[0,1/100],\quad\zeta\in[11,13],       \tag{5}
\]

with `1-T_zeta>0`, positive lapse, `F<0`, and negative **total** `b_v`.
In the formulas below, partial derivatives first treat `v,zeta` as
independent, with `c` fixed; total derivatives are then taken on the branch.
The denominator `D` vanishes at the corner of the extended box, but it is
not divided by in (4). All functions used to obtain the branch jets are
analytic through that corner. The literal coefficients `b1=Q_profile/D`
are still undefined there.

## 3. Exact interval Taylor algebra

`jets.py` implements the ring of two-variable polynomials truncated at
total degree three. A coefficient with index `(i,j)` encloses the actual
partial derivative divided by `i!j!` at every center in (5). Products use
the finite coefficient convolution and discard only terms with degree
above three. These are identities of jets, not a finite-increment
approximation to a function value.

For an analytic function jet `a`, let `a0` be its nonzero constant and write
formally `a=a0(1+h)`, where `h` has **exactly zero constant coefficient**.
In the truncated ring `h^4=0`, so

\[
a^p=a_0^p\sum_{n=0}^3\binom pn h^n.                                \tag{6}
\]

This proves the reciprocal and cube-root derivative recurrences without
a convergence assumption. In interval implementation, the constant is
removed by its coefficient index, not by interval subtraction of two
equal uncertain constants. Otherwise a spurious nonzero constant interval
would invalidate this nilpotent argument.

Every endpoint operation uses exact integer/Fraction interval arithmetic.
The reciprocal requires a constant interval excluding zero. The positive
cube root requires an explicitly proved rational bracket: its lower cube
must not exceed the entire base interval's lower endpoint, and its upper
cube must not be smaller than the upper endpoint. The S6.33 bracket
`199/100<b<201/100` passes this test on (5). At an exact center with `b^3=8`,
the zero-width bracket `[2,2]` is also checked directly. No numerical root,
floating point, finite difference, or sampled-grid inference is used.

Differentiating the explicit rational expressions for `L,H` in (4) uses
the necessary higher derivatives of `J,D,R` automatically. It does not
borrow an unspecified bound on their derivatives. The Taylor algebra is
of degree three in **these explicit expressions**, which already contain
the indicated first derivatives of the primitive profiles.

## 4. Implicit differentiation and total derivatives

Write `z_j=d^j zeta/dv^j` and `B=1-T_zeta`. At the fixed point,

\[
\begin{aligned}
z_1&=\mathcal T_v/B,\\
z_2&=(\mathcal T_{vv}+2\mathcal T_{v\zeta}z_1
                         +\mathcal T_{\zeta\zeta}z_1^2)/B,\\
z_3&=[\mathcal T_{vvv}+3\mathcal T_{vv\zeta}z_1
       +3\mathcal T_{v\zeta\zeta}z_1^2+\mathcal T_{\zeta\zeta\zeta}z_1^3
       +3(\mathcal T_{v\zeta}+\mathcal T_{\zeta\zeta}z_1)z_2]/B.
\end{aligned}                                                        \tag{7}
\]

For any function `a(v,zeta)`, its first three total derivatives are

\[
\begin{aligned}
a'&=a_v+a_\zeta z_1,\\
a''&=a_{vv}+2a_{v\zeta}z_1+a_{\zeta\zeta}z_1^2+a_\zeta z_2,\\
a'''&=a_{vvv}+3a_{vv\zeta}z_1+3a_{v\zeta\zeta}z_1^2
       +a_{\zeta\zeta\zeta}z_1^3
       +3(a_{v\zeta}+a_{\zeta\zeta}z_1)z_2+a_\zeta z_3.
\end{aligned}                                                        \tag{8}
\]

The implementation applies (8) to `b` through third order and `F` through
second order. Let primes in the rest of this section mean these total
derivatives. Since `N=2b'/F`, one can calculate the required quantities
without hiding a lapse derivative:

\[
\begin{aligned}
A&=4Q_{\rm profile}b'/(Fb^2),\\
A'&=A\left[Q_{\rm profile}'/Q_{\rm profile}+b''/b'-F'/F-2b'/b\right],\\
\ell&=(\log k)'=3b'/b-b''/b'+F'/F,\\
\ell'&=3[b''/b-(b'/b)^2]-[b'''/b'-(b''/b')^2]
                              +[F''/F-(F'/F)^2],\\
\frac{(\sqrt{k})_{uu}}{\sqrt{k}}&=\ell+2v\ell'+v\ell^2.
\end{aligned}                                                        \tag{9}
\]

The last identity uses `v=u^2`; replacing a partial derivative by a total
one only after this step would give the wrong operator. Every denominator
in (7)-(9) is explicitly enclosed away from zero. A whole-box Fraction
evaluation gives, with directed outward rounding,

| Quantity | Lower endpoint | Upper endpoint |
| --- | --- | --- |
| `A_v` (total) | `-128747179/1000000` | `37954017/1000000` |
| canonical pump | `-21449173/1000000` | `-2841831/250000` |

Thus `|A_v|<129` and `|pump|<22` on all of (5), not just on sampled fixed
points. The larger box includes the actual fixed-point graph. The verifier
can also replay all implicit and total intermediate jet enclosures.

## 5. Uniform pole cancellation

Set `Z=delta+8v`. Polynomial arithmetic gives

\[
D=Z-v^2E(v),\qquad
E(v)=\frac{20+40v+30v^2+8v^3}{d^4}.                                 \tag{10}
\]

All numerator coefficients are positive. Moreover,

\[
20d^4-(20+40v+30v^2+8v^3)=40v+90v^2+72v^3+20v^4\ge0,
\]

so `0<E<=20` for `v>=0`. The other useful inequality follows without any
small-parameter expansion:

\[
d^4D-Z=\delta(d^4-1)+12v^2+8v^3+2v^4\ge0,
\qquad D\ge Z/d^4.                                                   \tag{11}
\]

At `v=0`, the exact branch gives

\[
b=2,\quad b'=-c(c+2),\quad F=-4(c+2)/c,\quad Q_{\rm profile}=32/c.
\]

Consequently `A(0,c)=8c=16+8delta` exactly. Integrating the already-proved
bound on the **total** `A_v` along the actual branch at fixed `delta` gives

\[
|A(v,\delta)-16|\le8\delta+129v.
\]

The segment of integration lies entirely inside the proved branch box.
No bound on an inverse clock or a `delta` derivative is needed. Therefore

\[
\begin{aligned}
\left|\frac A D-\frac{16}Z\right|
&\le\frac{|A-16|}D+\frac{16v^2E}{DZ}\\
&\le d^4\max(8,129/8)+5d^4.                                        \tag{12}
\end{aligned}
\]

Here the second term uses `Z>=8v`; it is zero at `v=0`, where `delta>0`
keeps every displayed literal denominator nonzero. Finally,

\[
d^4\le(10001/10000)^4<1001/1000,
\]

and (9), (12) imply

\[
\left|V_\delta-16/Z\right|
<\frac{1001}{1000}\left(\frac{129}{8}+5\right)+22
=\frac{345169}{8000}<44.                                            \tag{13}
\]

The final strict rational margin is `6831/8000`. Every polynomial and
rounding direction needed for this calculation is separately checked.

## 6. Corner controls: bounded does not mean jointly continuous

At the excluded corner `(v,delta,zeta)=(0,0,12)`, exact rational jets of
the extended branch give

\[
\begin{gathered}
\zeta'=204,\quad\zeta''=8808,\quad\zeta'''=670824,\\
b'= -8,\quad b''=-56,\quad b'''=-4800,\\
F=-8,\quad F'=-24,\quad F''=-2080,\quad
Q_{\rm profile}=16,\quad Q_{{\rm profile},v}=-240.
\end{gathered}
\]

In particular, `A=16`, `A_v=-48`, `ell=-16`, and `ell_v=-432`.
The direct reconstruction `k=b^3/(2b'/F)` independently verifies the log
derivatives and the pump. Joint analyticity of these **regular** quantities
and the exact center identity imply

\[
A=16+8\delta-48v+O((v+\delta)^2),\quad
D=\delta+8v-20v^2+O(v^3),\quad
\text{pump}=-16+O(v+\delta).
\]

For finite `delta/v -> lambda>=0`, the pole-subtracted potential therefore
has the directional limit

\[
24-\frac{112}{\lambda+8}+\frac{320}{(\lambda+8)^2}.                  \tag{14}
\]

It tends to **15** along the admitted path `delta=v^2`, and to **24** along
`v=0, delta->0+`. The error terms above vanish after the indicated divisions
because `D` is bounded below by a positive multiple of `delta+v`.

Thus the combined remainder has no joint continuous extension at the
corner. Equation (3) must not be promoted to a uniform derivative bound for
that combined remainder. The corner calculation does not create a literal
`delta=0` action, and the two approach paths can both be taken with
strictly positive literal `delta`.

## 7. Comparison equation and remaining boundary

For the pole-only homogeneous equation, set `x=u/sqrt(delta)`:

\[
y_{xx}+\frac{16}{1+8x^2}y=0.
\]

The exact substitution `sqrt(8)*x=sinh(t)`, `y=sqrt(cosh(t))*psi` gives

\[
\psi_{tt}+\left[\frac74+\frac34\operatorname{sech}^2t\right]\psi=0,
\]

with log-wave frequency `sqrt(7)/2`. This substitution is an operator
identity, not a source or state choice. To turn the finite bound (3) into
a physical transfer result one must still specify the canonical endpoint
frames, initial/boundary data, the weighted source in (2), and a quantitative
propagation comparison. None is supplied merely by an instantaneous mass
or by the absence of a divergent pole-subtracted coefficient.

All new tests and computations here are local to this child. The frozen
action, source couplings, clock, and S6.33 branch remain unchanged. No
full-parent, fixed-low-frequency EFT, or UV conclusion is claimed.
