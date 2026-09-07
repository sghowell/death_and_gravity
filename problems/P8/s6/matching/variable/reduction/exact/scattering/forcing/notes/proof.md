# S6.36: a fixed prescribed pulse and separated actual own-f readouts

The input is one fixed smooth physical-metric tensor q, with **zero hidden
initial data**. It is not a matter source and is not prepared hidden data.
The result concerns the same off-shell `g=eta`, `theta=T` probe and the same
positive exact own-f background as S6.33. It compares the stated family
`c=2+delta` as delta tends to zero; the corner is not declared to be a
literal parent action. Constant M and tau are held fixed.

The pinned S6.35 report is
`c77d6d3e41d81c5d2bda33aeba201fe438620d872f3e9f2d49fde5469194641e`.
Its source normalization, canonical remainder estimate, and coefficient-16
connection are used. No coefficient-80 response result is imported.

## 1. Literal source, domain, and zero data

Set

\[
u=T/\tau,\quad L=1/100,\quad a=L/2,\quad0<\delta\le10^{-6}.
\]

In the norm-two tensor convention the physical quadratic action, with
the common positive factor \(M^2/(4\tau)\), is

\[
\int du\,[q_u^2+kQ_u^2-2b_1b(Q-q)^2],\qquad
k=b^3/N>0,\quad b_1=\tau^2\beta_1/M^2.
\tag{1}
\]

The actual own-f equation is
\((kQ_u)_u+2b_1b(Q-q)=0\). Define

\[
Y=\sqrt{k}\,Q,\qquad \Omega^2=2b_1b/k=\mathcal A/D,
\quad\mathcal A=2Q_{\rm profile}N/b^2.
\]

It becomes, without dropping either the pump or source weight,

\[
Y_{uu}+V_\delta Y=\sqrt{k}\,\Omega^2q,
\qquad V_\delta=\Omega^2-(\sqrt{k})_{uu}/\sqrt{k}.
\tag{2}
\]

The retained source is \(\sqrt{k}\Omega^2q\), not \(\sqrt{k}V_\delta q\).
The initial conditions are \(Q(-L)=Q_u(-L)=0\), equivalently
\(Y(-L)=Y_u(-L)=0\). The derivative in these conditions is
\(Q_u=\tau Q_T\).

S6.33 and S6.35 give, on the whole required domain,

\[
19/10<b,N<21/10,\quad
\left|V_\delta-\frac{16}{\delta+8u^2}\right|<44,\quad
\left|(\sqrt{k})_{uu}/\sqrt{k}\right|<22.
\tag{3}
\]

On the fixed loading interval \([-L,-a]\), these imply

\[
1<\frac{6859}{2100}<k<\frac{9261}{1900}<5,
\qquad0<V_\delta<\frac8{L^2}+44=80044,
\tag{4}
\]

and

\[
\Omega^2>
\frac{16}{10^{-6}+8L^2}-44-22
=\frac{15947134}{801}>\frac1{L^2}.
\tag{5}
\]

Thus the exact canonical source in (2) obeys
\(\sqrt{k}\Omega^2q\ge q/L^2\) for nonnegative q. This is derived
from the physical source map, not an assumed controllability coefficient.

## 2. One explicit delta-independent smooth pulse

Let \(S(x)\) be zero for \(x\le0\), one for \(x\ge1\), and

\[
S(x)=\frac{e^{-1/x}}{e^{-1/x}+e^{-1/(1-x)}}\quad(0<x<1).
\]

For one fixed \(\eta>0\), set

\[
q(u)=\eta S\!\left(\frac{u+95L/128}{L/64}\right)
             S\!\left(\frac{-81L/128-u}{L/64}\right).
\tag{6}
\]

Its support closure lies strictly inside \((-3L/4,-5L/8)\), with margin
\(L/128\) at each side. It is between zero and eta, and its plateau
has length \(5L/64>L/16\). Consequently

\[
\int q\,du\ge5\eta L/64>\eta L/16.
\tag{7}
\]

No transcendental quadrature is used for this area bound. Smoothness at
the ends follows from
\(d^n e^{-1/x}/dx^n=e^{-1/x}P_n(1/x)\), with polynomial recurrence
\(P_{n+1}(z)=z^2(P_n-P_n')\); every fixed derivative tends to zero.
The denominators in the interior definition of S are positive. In
particular q and all its derivatives vanish near both loading endpoints
and throughout the subsequent center-crossing interval. Formula (6)
contains no delta. Its physical duration and shape are fixed in T, not
in a shrinking inner coordinate.

The rest of the proof applies to any one fixed smooth \(0\le q\le\eta\)
supported inside \((-3L/4,-5L/8)\) with \(\int q\,du\ge\eta L/16\);
the displayed pulse proves that this generic fixed-input class is nonempty.

## 3. Positive loading proved by a first-zero argument

Let \(G(u,s)\) solve \(G_{uu}+V_\delta G=0\), with
\(G(s,s)=0\), \(G_u(s,s)=1\). For source-to-load lags
\(h=u-s\le L/4\), a first-zero bootstrap on \(G_u\) gives

\[
0\le G\le h,\qquad
G_u=1-\int_s^u V_\delta(v)G(v,s)dv
\ge1-80044h^2/2.
\]

At a putative first zero of \(G_u\), the last expression is at least
\(59989/80000>0\), a contradiction. Integrating the same lower bound
therefore proves throughout this lag interval

\[
G_u\ge\frac{59989}{80000}>\frac23,
\qquad
\frac G h\ge\frac{219989}{240000}>\frac9{10}.
\tag{8}
\]

These are uniform finite-delta bounds. For each point in the support (6),
the lag to \(-a\) is strictly between \(L/8\) and \(L/4\). The actual
zero-data Duhamel response and (5)--(8) yield

\[
Y_\delta(-a)\ge\frac{219989}{30720000}\eta
 >\frac9{1280}\eta,
\qquad
Y_{\delta,u}(-a)>\frac{\eta}{24L}.
\tag{9}
\]

The first displayed stronger constant already uses only the relaxed
area \(\eta L/16\) and lag \(L/8\). Keeping that stronger uniform
constant is important: a merely pointwise strict inequality would not
justify a strict inequality after taking a limit.

## 4. Convergence is needed only away from the center

On \([-L,-a]\), the S6.33 analytic joint background and all required
derivatives extend to delta zero. The potentially singular literal
denominator has the explicit margin

\[
D=2+\delta-\frac2{(1+u^2)^4}
\ge2-\frac2{(1+a^2)^4}>0.
\tag{10}
\]

The kinetic coefficient is also uniformly positive. Hence the actual
coefficients and source in (2) converge uniformly on this fixed compact
interval to their punctured delta-zero limits. Writing (2) as a first-order
system, the integral equation and Gronwall inequality prove convergence of
its zero-data solution there. This is ordinary compact-interval parameter
dependence. It does not assume a limit of the combined canonical remainder
at \((u,\delta)=(0,0)\).

For the center crossing put

\[
r_\delta=\sqrt{u^2+\delta/8},\quad
t=\operatorname{arsinh}(\sqrt8u/\sqrt\delta),\quad
Y=\sqrt{r_\delta}\,\psi,\quad Z=(\psi,\psi_t/\rho)^T,
\quad\rho=\sqrt7/2.
\]

At \(-a\), let \(w_\delta=Z(-a)\). Its limit is

\[
w=\begin{pmatrix}
Y_0(-a)/\sqrt a\\
[\sqrt a\,Y_{0,u}(-a)+Y_0(-a)/(2\sqrt a)]/\rho
\end{pmatrix},\qquad
\|w\|\ge w_1>\frac{9\eta}{1280\sqrt a}>0.
\tag{11}
\]

The strict inequality follows from the stronger uniform constant in (9).
Thus no hidden-data or source-loading premise has been added.

## 5. Recomputed finite-width transfer estimate

Between \(-a\) and \(+a\), q vanishes. The exact equation is the S6.35
homogeneous equation with its full bounded canonical remainder. Let
\(t_*=\operatorname{arsinh}(\sqrt8a/\sqrt\delta)\). Repeating S6.35's
finite-window comparison at this **smaller** half-width gives

\[
\int_{-t_*}^{t_*}|r_\delta^2R_\delta|dt
 <44(a^2+2a/2000)=\frac{33}{25000},
\qquad
\int_{|t|>t_*}\frac34\operatorname{sech}^2t\,dt
 \le\frac{3\delta}{32a^2}\le\frac3{800}.
\]

Using \(\rho>5/4\), the full perturbation budget is therefore

\[
J<\frac45\left(\frac{33}{25000}+\frac3{800}\right)
 =\frac{507}{125000}<\frac1{200}.
\tag{12}
\]

The reference flow bound \(B_0<2\), its one-\(B_0\) Duhamel estimate,
and the exact unitary endpoint wave frames give

\[
\|F_\delta-F_{\rm ref}(\theta_\delta)\|
 <e_*=\frac4{399},\qquad
\theta_\delta=\rho\operatorname{arsinh}(\sqrt8a/\sqrt\delta).
\tag{13}
\]

Here F is evolution in the actual Z variables. The free exterior used
in obtaining (13) is only a mathematical definition of finite-endpoint
coordinates. It is not a physical vacuum continuation.

The error in (13) is uniform, but is not assumed to tend to zero. No
limit, derivative bound, or outer asymptotic expansion of the central
remainder is used below.

## 6. Real phase matrix and two separated neighborhoods

Let the S6.35 Jost coefficients be \(A_J\) and \(B_J=i\beta\), where
\(\beta=\operatorname{csch}(\pi\rho)>0\) and
\(|A_J|^2=1+\beta^2>1\). In the declared plus/minus time-frequency
convention,

\[
T_t=\begin{pmatrix}\bar A_J&-\bar B_J\\-B_J&A_J\end{pmatrix},
\quad
V_0=\frac1{\sqrt2}\begin{pmatrix}1&1\\i&-i\end{pmatrix},
\quad D(\theta)=\operatorname{diag}(e^{i\theta},e^{-i\theta}).
\]

Since the negative endpoint has phase \(-\theta\), the Cauchy reference
is \(F_{\rm ref}(\theta)=V_0D(\theta)T_tD(\theta)V_0^*\), with a
D on **both** sides, not a unitary conjugation. If
\(A_Je^{-2i\theta}=\alpha_\theta+i\gamma_\theta\), it is the real map

\[
F_{\rm ref}(\theta)=
\underbrace{\begin{pmatrix}0&-\beta\\-\beta&0\end{pmatrix}}_{C_B}
+\underbrace{\begin{pmatrix}\alpha_\theta&-\gamma_\theta\\
\gamma_\theta&\alpha_\theta\end{pmatrix}}_{R_A(\theta)}.
\tag{14}
\]

Increasing theta by \(\pi/2\) reverses \(R_A\) and leaves \(C_B\)
unchanged. Also \(R_A^TR_A=|A_J|^2I\), so the two reference outputs on
the one fixed nonzero w differ in norm by exactly \(2|A_J|\|w\|\).
This effect would remain if B were artificially zero; it is not inferred
merely from a nonzero reflection coefficient.

Choose any fixed \(\theta_0\ge0\) and define, for integers \(n\ge2\),

\[
\delta_n^+=\frac{8a^2}{\sinh^2[(\theta_0+\pi n)/\rho]},\qquad
\delta_n^-=\frac{8a^2}{\sinh^2[(\theta_0+\pi n+\pi/2)/\rho]}.
\tag{15}
\]

Both sequences are positive and tend to zero, and they give exactly the
two phases in (14). Their domain has a rational proof: \(\pi>3\),
\(\rho<4/3\), and \(n\ge2\) imply
\(\sinh[(\theta_0+\pi n)/\rho]>\sinh(9/2)>315/16\), by the first two
positive terms of the sinh series. Thus both deltas are below
\(8a^2(16/315)^2<10^{-6}\). One direct proof of \(\pi>3\) integrates
\(1/(1+x^2)=\sum_{j=0}^7(-x^2)^j+x^{16}/(1+x^2)\) on \([0,1]\), giving
\(\pi>135904/45045>3\).

Uniform boundedness of the reference and actual maps, together with
\(w_\delta\to w\), shows that the actual outgoing vectors lie within
\(e_*\|w\|+o(1)\) of the respective reference centers
\((C_B\pm R_A(\theta_0))w\). The two error terms need not converge.
Consequently

\[
\liminf_{n\to\infty}
\|Z_{\delta_n^+}(+a)-Z_{\delta_n^-}(+a)\|
\ge(2|A_J|-8/399)\|w\|
>\frac{790}{399}\|w\|.
\tag{16}
\]

This proves separated actual neighborhoods, not equality of actual
subsequential limits to the reference centers. The domain bound n>=2
does not assert that the asymptotic separation already holds at n=2.

## 7. A separated original-Q scalar readout

There is no need to stop at a chosen two-vector norm. Represent
\(w=(w_1,w_2)\) by \(w_1+iw_2\), and choose

\[
\theta_0=\tfrac12\arg[A_J(w_1+iw_2)]\pmod\pi,\qquad0\le\theta_0<\pi.
\tag{17}
\]

The product is nonzero by (11), so this is a fixed, delta-independent
choice. Then \(R_A(\theta_0)w=(|A_J|\|w\|,0)^T\). The first components
of the reference centers are separated by \(2|A_J|\|w\|\), and the same
two error bounds apply to that component.

The full actual physical-to-Z map at \(u=\sigma a\) is

\[
Z=\begin{pmatrix}r_\delta^{-1/2}&0\\
-u/(2\rho r_\delta^{3/2})&\sqrt{r_\delta}/\rho\end{pmatrix}
\begin{pmatrix}\sqrt{k}&0\\k_u/(2\sqrt{k})&\sqrt{k}\end{pmatrix}
\binom Q{Q_u}.
\tag{18}
\]

It is convergent and invertible at either fixed nonzero endpoint, with
determinant k/rho. All of \(k_u\) is retained in both the loaded state
and the two-vector map. Its first row gives the exact scalar identity
\(Q(+a)=\sqrt{r_\delta(+a)/k_\delta(+a)}\,Z_1(+a)\).

Let \(k_0=k_0(+a)>0\). The scalar prefactor tends to
\(\sqrt{a/k_0}\). The outgoing Z vectors are bounded, so replacing this
prefactor by its limit changes the readout only by \(o(1)\). Equations
(11), (16)--(18), and \(k_0\le9261/1900<5<25/4\) prove

\[
\begin{aligned}
\liminf_n[Q_{\delta_n^+}(+a)-Q_{\delta_n^-}(+a)]
&\ge(2|A_J|-8/399)\sqrt{a/k_0}\,\|w\|\\
&>\frac{790}{399}\frac9{1280}\frac25\eta
=\frac{237}{42560}\eta>\frac{\eta}{200}.
\end{aligned}
\tag{19}
\]

In particular the bounded family has

\[
\limsup_{\delta\downarrow0}Q_\delta(+a)
-\liminf_{\delta\downarrow0}Q_\delta(+a)>\eta/200.
\tag{20}
\]

Thus neither the scalar Q readout nor the actual pair \((Q,Q_u)\) is
Cauchy as delta tends to zero. No assertion about the separate derivative
readout's limsup/liminf is needed.

## 8. Scope and relevance to matching

This is one fixed smooth prescribed-q experiment, with a proved nonzero
loading from zero hidden data, on a fixed physical interval from
\(T=-\tau/100\) to \(T=\tau/200\). It strengthens a prepared-data or
shrinking-pulse result: neither the input nor its support depends on delta.
It does not assume convergence of the central remainder or promote a
reference phase curve to an exact endpoint limit.

Everything is at the linear tensor level about the specified off-shell
probe. Maintaining q may require driving the physical metric. It is not
an on-shell g solution, a free-chi or conserved matter-source response, a
fixed physical-frequency EFT band, a nonlinear error estimate, or a
vacuum/UV statement. Taking eta small scales this linear result; it does
not by itself prove uniform nonlinear validity as delta tends to zero.

For the adopted B matching gate this supplies a source-preserving warning
against claiming a delta-independent, memory-free own-f response on this
particular input class. A matching proposal that restricts its frequency
domain or prepares another state needs its own quantitative comparison.
The result supplies neither the adopted V vacuum test nor the G dispersive
remainder, and it does not exclude all EFT reductions or close original
P8. Further refinements of this same probe are not a substitute for an
actual common-parent matter-frame matching construction.
