# S6.26.FORCED — causal nonlocal physical response

This gate concerns the unchanged S6.20 VARIABLE local background and its full
linear two-tensor action, with physical matter/source metric `g`. It is a
causal **nonlocal** approximation theorem, not a local light-field EFT or an
original C/D matching verdict. The written proof, exact identities and
source-aware independent audit are replayed by a read-only certificate.

## Domain and state/source data

Use `u=T/tau`, `M,tau>0`, `c=2+delta`, `d=1+u^2`, `a=d^2`. For each fixed
`0<r<=1/100`, take `I=[-r,r]`, `0<delta<=min(10^-9,r^2/100)` and the whole
spatial band `0<=K=(tau*k_com)^2<=4`. All four physical tensor Cauchy data
vanish at `u=-r`. This is a specified classical linear-response preparation,
not selection of a quantum vacuum or of S6.23's analytic prepared sector.

A real external spatial TT stress has the source amplitude
`sigma=tau^2*Pi/M^2` from `+1/2 integral dT a^3 Pi gamma_g`. Its time profile
may be any `L-infinity(I)` function, with `S=||sigma||infinity`. Smooth profiles
flat in a left neighborhood are included. A single real Fourier polarization
is transverse and traceless; the source is conserved on the background. This
does not make it nonlinear matter or assert its backreaction is negligible at
arbitrary amplitude. The theorem is a linear transfer-operator estimate.

## Precise claim

Let `l,Q` be the normalized canonical fields. Define `L_L=d_u^2+A`,
`L_H=d_u^2+B`, `Bop=E+u*fc*d_u`, `Bopstar=C+u*dc*d_u`, with the exact literal
coefficients and source weights in `operator.py`. Both inverses below are
retarded with zero data at `-r`. The full equations are

`L_L l+Bopstar Q=jL*sigma`, `L_H Q+Bop l=jH*sigma`.

The exact elimination is the causal Schur equation

`(L_L-Bopstar*GH*Bop)l=jL*sigma-Bopstar*GH(jH*sigma)`.

Define `l0=GL(jL*sigma)` and
`Qapp=G0[jH*sigma-Bop*l0]`, where `G0` is the explicitly normalized retarded
inverse of `d_u^2+80/(delta+8u^2)`, with its exact finite-delta clock and memory.
The physical output is `gapp=ag*l0+bg*Qapp`, using the **actual** delta-dependent
physical map, not its limiting or locked version. Then

`||gamma_g-gapp||infinity <=1000*r^4*S`.

At `r=1/100` this is `10^-5*S` for every admitted delta. It is an absolute
source-to-field norm bound; neither a relative error where the response can
vanish, nor an error tending to zero with delta at fixed r, is asserted.

## Separate interpretation controls

The exact symplectic projection of the forced response onto S6.23's analytic
two-dimensional sector has source factor `a^3/(2*Keff)`, not `2`. Its retarded
diagonal derivative tends to `2/5` at the center as delta tends to zero, while
the actual full `g` diagonal derivative is `2`. This is a source/state control,
not a low-frequency omission-error lower bound.

For a nonzero `sigma in H1_0(I)`, define its RMS derivative temporal scale
`Omega_sigma=||sigma'||2/||sigma||2` and the instantaneous heavy-diagonal
stiffness proxy `m_proxy=inf_I sqrt(B)`. Both are explicitly specified proxies,
not a Fourier bandlimit or a proved rolling frequency gap/cutoff. The bound
`Omega_sigma/m_proxy>6/13` excludes a parametrically small version of this
particular derivative hierarchy as delta tends to zero. It does not exclude
other source classes, longer preparations, or all possible low-energy EFTs.

No characteristic, causality violation, positivity, UV, scalar/vector health,
global completion, nonlinear source, quantum particle-production, or original
S6/P8 closure claim follows. The memory kernel has not been replaced by a
local derivative expansion. All input sources and ancestor certificates stay
immutable; a future certificate must pin S6.23 and its S6.21/S6.20 lineage.
