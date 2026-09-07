# S6.21 RESPONSE: fixed-slice physical two-tensor transfer

This is a descendant of frozen `P8-S6.20.VARIABLE`, not a modification of its
action, physical metric, scalar sources or domain. The adopted S6 matching
contract is unchanged. The result is a classical linear-response theorem for
the **full coupled two-TT operator**, with its actual physical source/data map.
No scalar/vector health, low-temporal-band reduction or UV verdict follows.

## Fixed model, parameters and norms

Keep the VARIABLE action and its actual sourced local solution, with
`u=T/tau`, `M,tau>0`, `c=2+delta`, `a=(1+u^2)^2`, `b=2/(1+u^2)^2`.
The physical matter metric is `g`, `N_g=1`, `N_f=c`. Restrict

\[
0<\delta\le10^{-9},\quad 1\le\bar k=\tau k_{\rm com}\le2,
\qquad a_*:=1/100.
\]

Both observation slices `u=+-a_*` and the optional source interval
`[-2a_*,-a_*]` lie inside the parent's actual `|u|<=1/10` chart. Here `a_*`
is a slice radius, not the background scale factor `a(u)`. No global or
larger-duration extension is assumed. `delta=0` is used only to define a
punctured limiting operator, never as a regular parent through `u=0`.

Canonical variables below are normalized by `M`; derivatives are `d/du`.
Use Euclidean vector norms on real four-component data (or their complex
extension) and induced operator norms. A source uses the explicitly named
dimensionless `L1(du)` norm. Cauchy derivatives are not canonical momenta.

## Exact maps and matched theorem

Let `Z=(l,l',Q,Q')`. The exact physical map `C_delta(u)` sends it to
`G=(gamma_g,gamma_g',gamma_f,gamma_f')`; it includes derivatives of every
normalization and weight. These are actual amplitudes and `tau` times their
proper-time derivatives. For raw proper derivatives multiply rows 2 and 4
by `tau^-1`. Define, on `u=sigma r`, `r>0`,

\[
Y=(l,l',Q/\sqrt r,(uQ'-Q/2)/(\mu\sqrt r)),\qquad \mu=\sqrt{39}/2.
\]

Write `Z=W_sigma(r)Y`. The **full** delta-zero punctured equation has
`Y_{log r}=(A_0+R_sigma(r))Y`, where `A_0` is zero on the first two
components and `[[0,mu],[-mu,0]]` on the last two. Its canonical outer wave
matrix is the unique Volterra solution

\[
\Psi_\sigma(r)=r^{A_0}\left[I+
 \int_0^r s^{-A_0}R_\sigma(s)\Psi_\sigma(s)\,ds/s\right].
\tag{1}
\]

All couplings and the fixed momentum remain in this matrix. The integral
converges absolutely since `||R_sigma(r)||<=64r`. Its series has a certified
factorial tail. Define the actual endpoint maps

\[
P_\delta^\sigma=C_\delta(\sigma a_*)W_\sigma(a_*)\Psi_\sigma(a_*).
\]

If `T_delta^phys` is the exact physical four-Cauchy-data propagator between
the fixed slices, then

\[
\left\|(P_\delta^+)^{-1}T_\delta^{\rm phys}P_\delta^-
       -\operatorname{diag}(I_2,S_{\sqrt\delta}^{\rm real})\right\|
\le40000\delta^{1/3}.
\tag{2}
\]

The endpoint maps are uniformly invertible: `||P_delta^sigma||<3000` and
`||(P_delta^sigma)^-1||<280` in this dimensionless physical clock. Thus (2)
is not a quotient by an unobservable field. Exact maps can instead be used
for a less enlarged physical norm. The error is at most `1/250` when
`delta<=10^-21`; it is **40**, not small, at the upper box endpoint `10^-9`.

Here `S_epsilon` in the complex normalized power basis
`|u|^(1/2+-i mu)` is

\[
S_\epsilon=
\begin{pmatrix}
-\overline B&\overline A(4\sqrt2/\epsilon)^{2i\mu}\\
A(\epsilon/(4\sqrt2))^{2i\mu}&-B
\end{pmatrix},\quad
A=\frac{\Gamma(1-i\mu)\Gamma(-i\mu)}
        {\Gamma(3/2-i\mu)\Gamma(-1/2-i\mu)},\quad
B=\frac{i}{\sinh(\pi\mu)}.
\tag{3}
\]

The normalized power basis is unitarily equivalent to the real cosine/sine
basis of (1). In particular `|A|^2-|B|^2=1`, `0<|B|<1/1000`, `||S||<2`.
Its determinant is `-1` because radial orientation reverses, not a failure
of physical Cauchy symplecticity. Even setting `B=0` leaves the transmitted
`A` phase; nonconvergence is not attributed solely to reflection.

## Physical source and prepared-state controls

For a unit-norm transverse polarization and fixed nonzero `kbar`, the
linear external TT stress `delta T^{ij}=a^-2 Pi(u)e^{ij}cos(k_com z)`,
`delta T^{0mu}=0`, is covariantly conserved for arbitrary smooth `Pi`.
The probe action is `1/2 integral dT a^3 Pi gamma_g`. Put
`sigma=tau^2 Pi/M^2`. Its canonical forcing is

\[
j_l=\frac{a^3\sigma}{2f_\Sigma},\qquad
j_Q=-\frac{a^3w_2\sigma}{2f_R};\qquad
(j_l,j_Q)_{u=0,c=2}=(1,-2)\sigma/\sqrt5.
\]

For every fixed smooth compact source on `(-2a_*,-a_*)`, with zero retarded
data on its left, its incoming outer coefficient differs from the
delta-zero coefficient by at most

\[
3\,000\,000\delta\,\|\sigma\|_{L^1(du)}.
\tag{4}
\]

An explicit `g`-only loading operator constructs such a source for either
unit target `C=(0,0,1,0)` (incoming relative mode) or `C=(1,0,0,0)` (prepared
regular-light mode); the target is a full homogeneous punctured solution,
not a locked ansatz. The same fixed source is then used at every positive
delta. Its exact finite norm `S=integral|sigma|du` remains visible. The
outgoing error against `diag(I,S_epsilon)C` is bounded by

\[
40000\delta^{1/3}\|C\|
+(2+40000\delta^{1/3})3\,000\,000\delta S.
\tag{5}
\]

Thus prepared-light leakage has an explicit vanishing upper bound, not a
claimed sharp `O(delta)` rate. A prepared nonzero relative component gives
phase-dependent fixed-slice physical response. After the source ends the
map to the **g-only four-jet** `(g,g',g'',g''')` is invertible with determinant
`(U_TT/K_g)^2>0`. We do not assume the two-component projection `(g,g')`
alone cannot cancel, and a finite jet is not a low temporal-band observable.

## Scope

The source is an infinitesimal external probe; no nonlinear source energy
condition, tensor backreaction or vacuum-production conclusion is imposed.
The time support is fixed, but a compact smooth pulse is not band-limited.
The fixed `kbar` spatial band is not a low **temporal** band. No source/state
selection is supplied by the adopted low-energy matching contract here.
In particular the estimates do not compare a proven omission error below a
cutoff to the center's `O(delta)` locked-cone excess. No uniform rolling gap,
full light-only EFT exclusion or original C/D operator match is certified.
Original P8(b), S6 and P8 remain open.
