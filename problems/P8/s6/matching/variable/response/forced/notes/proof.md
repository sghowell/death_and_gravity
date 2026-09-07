# Actual-source causal memory reduction

This proof is for the reviewed S6.26 claim in `FORMULATION.md`. It does not
certify a local EFT by naming the instantaneous diagonal coefficient a mass.
Every estimate is uniform on the stated finite positive-delta domain.

## 1. Full physical operator, source, and initial data

Write `u=T/tau`, `c=2+delta`, `d=1+u^2`, `W=c*d^12+8`, `a=d^2`, and
`K=(tau*k_com)^2`. Primes in this document are `d/du` unless a different
coordinate is specified. The canonical fields of the pinned parent are
divided by `M`. With

```
theta = 6u(c d^12-8)/(d W),
omega = -24 sqrt(2c) u d^5/W,
w1 = c d^12/W,  w2 = 8/W,
cf2 = c^2 d^8/4,  sw = sqrt(8c) d^6/W,
cross = (K/d^4) sw (cf2-1),
D = c d^4-2,
N = 16(1-u^2)[8/(c d^14)+1/d^2],
```

the full equations, including the moving-weight terms, are

```
l''+A l+C Q+u dc Q' = jL sigma,
Q''+B Q+E l+u fc l' = jH sigma,                         (1)
A = (K/d^4)(w1+cf2 w2)-theta'-theta^2,
C = cross-2omega'-2omega theta, dc=-2omega/u,
E = cross-2omega theta,         fc= 2omega/u,
B = N/D+b_an,
b_an = (K/d^4)(w2+cf2 w1)+theta'-theta^2-4omega^2.
```

The expressions `dc,fc` have their regular values at zero. The multiplication
and derivative parts are related by
`dc=-fc`, `C=E-(u fc)'`; hence `Bopstar` is the formal adjoint of
`Bop=E+u fc d_u`, with the boundary term `u fc*l*Q`. Neither is omitted.

The actual physical map and source are

```
g=ag l+bg Q,
ag=2 sqrt(c) d^3/sqrt(W),  bg=-4 sqrt(2)/(d^3 sqrt(W)),
jL=sqrt(c) d^9/sqrt(W),   jH=-2 sqrt(2) d^3/sqrt(W),
jL=a^3 ag/2,             jH=a^3 bg/2,
sigma=tau^2 Pi/M^2.                                      (2)
```

They follow by varying `+1/2 integral dT a^3 Pi*g`, not by choosing a source
in canonical heavy coordinates. The exact identity `ag*jL+bg*jH=2` is the
physical `g` retarded diagonal acceleration normalization.

For a real spatial plane wave choose a constant polarization `e_ij` with
`k_i e_ij=0`, `e_ii=0`, and stress `T_ij=a^2 Pi(T)e_ij cos(k.x)`, `T_0mu=0`.
On the homogeneous background its spatial divergence vanishes by
transversality and its temporal divergence is proportional to the zero trace.
Thus the external linear TT source is conserved, for arbitrary time profile.
`K=0` is also allowed. Normalization of the chosen polarization is the one in
the pinned TT action and (2); a rescaling rescales `Pi` and `sigma` together.

For each admitted positive delta all coefficients are smooth on the compact
interval. Zero physical Cauchy data at `-r` are equivalent to zero canonical
data by the invertible exact physical Cauchy map. An `L-infinity` source gives
a unique `W2,infinity` solution (and hence a `C1` solution); smooth sources
give smooth solutions. Finite spatial Fourier superpositions obey linearity,
but no unmentioned integration norm over momentum is claimed.

## 2. Continuous coefficient enclosures, including K=0

Here `0<r<=1/100`, `t=r^2`, `v=u^2<=t`, `0<delta<=t/100`, and `0<=K<=4`.
Set `cp=2+t/100`, `dp=1+t`. Then `2<=c<=cp`, `1<=d<=dp`, `W>=10`, and
`c d^12<8`. No lower bound on K enters this section.

Write `theta=u*T(v)` and `omega=u*O(v)`, differentiating `v` at fixed c.
Directly from the literal rational functions,

```
|T| <=18/5 = T0,
|O| <= (48+12t/100) dp^5/10 = O0,
|T_v| <=6[16*12cp dp^11/100+3/5] = T1,
|O_v| <=O0[5+12cp dp^11/10] = O1.                       (3)
```

For example `T=6(P-8)/(d(P+8))`, `P=c*d^12`, so
`T_v=6[16P'/(d(P+8)^2)-(P-8)/(d^2(P+8))]`. Also
`O_v=O[5/d-12c d^11/(P+8)]`. The square-root bound used in (3) is
`sqrt(2c)<=2+delta/2`, obtained by squaring positive sides.

Put

```
Mtheta=T0+2t T1+t T0^2,
SW=(4+t/100)dp^6/10,
Z=(cp^2 dp^8/4-1)/t,
X=4 SW Z.
```

The removable `Z` is a polynomial in t with nonnegative coefficients. The
following are upper bounds with nonnegative-power dependence on t:

```
|A| <=4(cp dp^12+2cp^2 dp^8)/10+Mtheta <9,
|C| <=t X+2(O0+2t O1)+2t O0 T0 <12,
|E| <=t(X+2O0 T0) <60r^2,
|dc|,|fc| <=2O0 <10,
|b_an| <=4(8+cp^3 dp^20/4)/10+Mtheta+4t O0^2 <9.       (4)
```

All strict inequalities are checked rationally at `t=1/10000`, which is
sufficient by their displayed positive coefficients. They are not numerical
samples of the actual operator. In particular the small E bound is freshly
proved on the shrinking r-variable box, not inferred from a K>=1 theorem.

For the pole isolate `D0=delta+8v`. One has exactly

```
D-D0=delta*v(4+6v+4v^2+v^3)+v^2(12+8v+2v^2)>=0,
D0^2=32delta*v+(delta-8v)^2.
```

On the larger rectangle `0<=v<=1/100`, `delta>=0`, the nonnegative N obeys
`N<=80`, `|N_delta|<=32`, `|N_v|<=1008`. Consequently
`|N-80|/D<=126`. Moreover the two coefficients in `D-D0` are bounded by
`41/10` and `121/10`. Thus

```
|N/D-80/D0|
 <=126+80[(41/10)/32+(121/10)/64]
 =1211/8 <152.                                            (5)
```

The derivatives in this argument are differentiated from the explicit N in
(1). Combining (4) and (5) gives

`B=80/(delta+8u^2)+b_delta`, `|b_delta|<161`.               (6)

Only `b_an=B-N/D` is analytic at the joint origin. The bounded `b_delta` in
(6) need not be jointly continuous there; no bound on its delta derivative
has been used.

The actual physical/source weights obey, by (2),

```
ag^2 <=4cp dp^6/10<1,       bg^2<=32/10<4,
jL^2 <=cp dp^18/10<1,       jH^2<=8dp^6/10<1.             (7)
```

This supplies `|ag|<1`, `|bg|<2`, and both source weights below one without a
limiting-field-map replacement.

## 3. Exact retarded pole kernel and uniformly bounded heavy memory

Let `GH` denote the zero-data inverse of `d_u^2+B`, and define

```
rho=sqrt(u^2+delta/8), z=asinh(sqrt(8)*u/sqrt(delta)),
Q=sqrt(rho)*psi, dz/du=1/rho, mu=sqrt(39)/2.
```

The exact heavy equation with forcing F is

`psi_zz+[mu^2+(3/4)sech^2 z+rho^2 b_delta]psi=rho^(3/2)F`. (8)

Here `3<mu<13/4`, `rho<=R=(101/100)r` on I. The upper R follows from
`rho^2<=r^2(1+1/800)<(101r/100)^2`. The canonical z-state is
`Z=(psi,psi_z/mu)`, with the ordinary Euclidean norm; induced operator norms
below use that same norm. The free mu generator is skew-adjoint. Factoring
out its rotations and applying the integral norm inequality gives, on every
ordered subinterval, propagator norm at most

```
exp[ integral |(3/4)sech^2 z+rho^2 b_delta|/mu dz ]
 <=exp(1/2+eta)<2,
eta=(322/3)r R <=(16261/150)r^2.                          (9)
```

We used `integral_R sech^2 z dz=2` and
`integral_I rho^2 dz=integral_I rho du<=2rR`. The same bound, with eta zero,
holds for the universal equation obtained by omitting only `rho^2 b_delta`.
The exact elementary upper bound
`exp(x)<=1+x+x^2/[2(1-x/3)]`, `0<=x<3`, follows from
`n!>=2*3^(n-2)` for `n>=2`; at the rational exponent in (9) it is below two.

For an explicit universal retarded kernel use

```
f+(z)=exp(i mu z) 2F1(-1/2,3/2;1-i mu;(1-tanh z)/2),
Gpsi0(z,w)=Im[f+(z) conjugate(f+(w))]/mu,  z>=w,
G0(u,s)=sqrt(rho(u)rho(s))*Gpsi0(z(u),z(s)), u>=s,         (10)
```

and set the kernel to zero for `u<s`. Its Jost normalization is solely a
basis choice for the universal auxiliary ODE, not an asymptotic state or an
assumption that the physical background extends beyond I. Substitution of
the hypergeometric differential equation gives the universal equation in (8). The defining
local series at its zero argument gives `f+(z)~exp(i mu z)` as `z->+infinity`.
The real-potential Wronskian is therefore
`W(f+,conjugate(f+))=-2i mu`. This proves the positive derivative jump in
(10). The coframe clock has not disappeared: `rho*z'=1` makes the Q
Wronskian unchanged and gives `partial_u G0(s+,s)=1`. For finite real z the
hypergeometric argument lies strictly in `(0,1)`; its c parameter is not an
integer pole. The single Jost representation uses no potentially degenerate
pair at hypergeometric infinity. [DLMF 15.10.1–2](https://dlmf.nist.gov/15.10)

For a bounded F the forcing in the z-state has magnitude
`rho^(3/2)|F|/mu`. Since `dz=du/rho`, (9) implies

```
||GH F||infinity <=(4/3)r R ||F||infinity,
||u (GH F)'||infinity <=(14/3)r R ||F||infinity.           (11)
```

For the second inequality use
`Q'=rho^(-1/2)[psi_z+(rho'/2)psi]`, `|rho'|<=1`, and
`|u|<=rho`. This proves the weighted derivative estimate without falsely
claiming a uniform unweighted Q' bound at the center. Equations (4),(11) give

```
||GH||Linfinity->C0 <=(101/75)r^2,
||Bopstar GH||Linfinity->Linfinity
 <=[12*(4/3)+10*(14/3)]rR <64r^2.                        (12)
```

Duhamel's identity between the full and universal z propagators has norm
at most `4eta`, from the two factors below two. A second retarded source
integration, with the same weights as (11), yields

```
||GH-G0||Linfinity->C0
 <=(8/3)eta*rR <=(16*161/9)r^2 R^2 <292r^4.             (13)
```

These are full finite-delta operator bounds. Neither the z interval's
increasing length nor a fixed-frequency WKB approximation is hidden in them.

## 4. Full causal elimination, light feedback and physical error

Let `GL` be the zero-data inverse of `d_u^2+A`. Define the Banach norm
`N(l)=||l||infinity+r||l'||infinity` on C1 functions with zero left data.
From the scalar Volterra series and `|A|<=9`, putting `L=2r`,

```
||GL F||infinity <=(L^2/2)cosh(3L)||F||infinity,
||(GL F)'||infinity <=L cosh(3L)||F||infinity.
```

The series comparison `(2n)!>=2^n n!` gives
`cosh(6r)<=exp(18r^2)<=1/(1-18r^2)`. Consequently

`N(GL F)<=4r^2 exp(18r^2)||F||infinity<5r^2||F||infinity`. (14)

Also `||Bop l||infinity<=(10+60r^2)N(l)` by (4). With no homogeneous
heavy solution in the specified initial data, exact elimination of (1) is

```
Q=GH(jH sigma-Bop l),
l=GL[jL sigma-Bopstar GH(jH sigma)+Bopstar GH Bop l].     (15)
```

The feedback operator in the second line has norm

`theta_feedback<=320r^4(10+60r^2)<1/10000`.               (16)

This is a uniformly convergent causal Neumann series, not a replacement of
the retarded inverse by `1/B`. It agrees with the unique original ODE
solution. From (7),(12),(14), with `S=||sigma||infinity`,

`N(l)<=5r^2(1+64r^2)S/(1-theta_feedback)<6r^2S`.         (17)

Define the approximation exactly by

```
l0=GL(jL sigma),
F0=jH sigma-Bop l0,
Qapp=G0 F0,
gapp=ag l0+bg Qapp.                                     (18)
```

The light Green function and all source/map coefficients in (18) are the
actual ones at delta; only the heavy bounded correction and feedback are
approximated. Since `N(l0)<=5r^2S`,

```
||F0||infinity <=(1+50r^2+300r^4)S <(101/100)S,
N(l-l0) <=320r^4(1+60r^2+360r^4)S <323r^4S.             (19)
```

Split the exact heavy error without suppressing the source term:

`Q-Qapp=(GH-G0)F0-GH Bop(l-l0)`.

It follows that

```
||Q-Qapp||infinity
 <={292*(101/100)+(101/75)r^2(10+60r^2)*323}r^4 S
 <296r^4S.                                               (20)
```

Every rational comparison in (9),(12)–(20) is monotone in r and is replayed
at `r=1/100`. Finally the actual physical map (7) gives

`||g-gapp||infinity<(323+2*296)r^4S=915r^4S<=1000r^4S`.   (21)

This is a fixed physical-duration interval `[-r*tau,r*tau]` when r,tau are
fixed. For `r=1/100`, its absolute error is at most `10^-5*S`. It does not
shrink with delta at fixed r, although the bound is uniform as delta tends
to zero. A shrinking r is a different limit requiring `delta<=r^2/100`.
This delta-independent upper bound alone cannot resolve an `O(delta)`
locked-cone effect or transfer the prepared center coefficient to the
arbitrarily forced problem; it does not prove the actual error stays large.
No relative-to-response lower denominator, derivative error, or positivity
of the response is claimed. Source norm and physical output normalization
cannot be dropped from this statement.

## 5. Source-free analytic preparation is not arbitrary-source closure

Let `z_e,z_o` be S6.23's exact prepared solutions, with canonical momenta
`pi_l=l'-2omega Q`, `pi_Q=Q'`, conserved symplectic pairing `Omega>0`, and
physical components `g_e,g_o`. At fixed K let
`Wg=g_e g_o'-g_e' g_o`, `Keff=Omega/Wg`. Their positive normalization is a
pinned source-free statement, not a choice of heavy retarded initial data.

For a full forced solution, project at each u onto the symplectic span of
these two homogeneous solutions. If the projected coefficients are
`alpha,beta`, the exact symplectic forcing identity and (2) give

```
alpha'=-a^3 g_o sigma/(2Omega),
beta'= a^3 g_e sigma/(2Omega).
```

For `gP=alpha*g_e+beta*g_o` the terms in its first derivative cancel. Its
second-order closed equation therefore has source

`gP''+Fprep*gP'+Gprep*gP = a^3 Wg sigma/(2Omega)`.

This is the correctly sourced **projection**, with coefficient
`a^3/(2Keff)`; it is not the full physical g equation. The prepared center
limit has `Omega->1`, `Wg->4/5`. Thus the projected retarded diagonal
derivative tends to `2/5`, whereas the exact full value in (2) is `2`. The
complementary value tends to `8/5`, not zero. This detects the unsupported
step of treating an arbitrary g-only source as if it stayed in the prepared
subspace. It does not assert that an impulse/coincident-time diagnostic
controls a temporally filtered low-frequency measurement.

In particular nonzero homogeneous initial data add their own causal response
to (15); they cannot silently be removed or identified with a quantum state.
An analytic retuned preparation may be a useful different input, but (21)
is for the stated zero-data forced problem only.

## 6. A limited RMS hierarchy obstruction, not a universal EFT no-go

For nonzero `sigma in H1_0(-r,r)` define

`Omega_sigma=||sigma'||L2/||sigma||L2`.

Zero extension has the corresponding Fourier second-moment interpretation.
This is a specified derivative/RMS temporal scale, not an exact bandlimit.
The interval Poincare inequality and `pi>3` give
`Omega_sigma>=pi/(2r)>3/(2r)`.

The instantaneous heavy-diagonal stiffness in (6) is positive throughout I:

```
r^2 B(u)>=8000/801-161r^2>9,
r^2 B(r)<=10+161r^2<(13/4)^2.
```

Consequently the *defined* proxy `m_proxy=inf_I sqrt(B)` obeys
`m_proxy<13/(4r)`, and

`Omega_sigma/m_proxy>6/13`.                              (22)

The ratio is unchanged by conversion to proper time, because both scales
acquire `1/tau`. Thus delta alone cannot produce a parametrically small
version of this declared derivative-versus-instantaneous-stiffness hierarchy
for these compact pulses on the fixed interval. This does not identify
`m_proxy` with a full rolling normal-mode gap, strong-coupling scale or
cutoff, and (22) is not a theorem excluding every low-energy source class.

As a separate elementary control, an L2 function with strictly compact
temporal Fourier support is entire by its finite-band Fourier integral;
if it vanishes on a past open interval, the identity theorem makes it zero.
Hence exact temporal bandlimiting and a nontrivial strictly causal pulse
cannot both be required. Approximate spectral concentration, persistent
forcing, longer preparation histories, and specially correlated initial
heavy data are distinct questions requiring new quantitative estimates.

None of these arguments supply an original C/D operator match, high-energy
positivity verdict, physical characteristic cone, or global UV completion.
