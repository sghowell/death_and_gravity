# P8: causal matter response and a global reciprocal boundary

This continues the [uniform scalar checkpoint](assessment-2026-09-07-p8-uniform-scalar-principal-sector.md).
Two further calculations separate what this parent actually proves from
what is still needed for original P8. Neither changes the original 32
linear row verdicts or the completed scoped photon objective. Original P8
remains open; no new user choice or permission is currently needed.

## A source-aware response, with its memory retained

The [forced-response gate](../problems/P8/s6/matching/variable/response/forced/FORMULATION.md)
keeps the same local action, physical g metric and time-dependent canonical
map. All four initial data vanish at u=-r. The external source is an
arbitrary bounded, background-conserved g-TT probe, with dimensionless
amplitude sigma=tau^2*Pi/M^2. The theorem covers

```
0<r<=1/100, |u|<=r,
0<delta=c-2<=min(10^-9,r^2/100), 0<=K=(tau*k_com)^2<=4.
```

There is no implicit choice of an analytic heavy state and no import of
the earlier transfer theorem outside its K>=1 band. The exact coupled
equations, written with their formal-adjoint derivative connections, are

```
L_L l + Bopstar Q = jL sigma,
L_H Q + Bop l    = jH sigma.
```

Both source components are fixed by the literal physical action. Their
combination gives the full retarded physical-g diagonal source factor 2.
Omitting the heavy source does not preserve that normalization.

The retarded heavy inverse G_H gives the exact causal Schur equation

```
(L_L-Bopstar G_H Bop) l = jL sigma-Bopstar G_H(jH sigma).
```

Its feedback is bounded by `320*r^4*(10+60*r^2)<10^-4` in the specified
light norm. An explicit approximation keeps the universal heavy kernel
G0 of `d_u^2+80/(delta+8u^2)`, with the actual finite clock
`z=asinh(sqrt(8)*u/sqrt(delta))`. The auxiliary hypergeometric solution
fixes a classical retarded Green function, not an asymptotic quantum
state or an extension of the physical background. The bounded remainder
in the heavy operator is controlled without falsely assuming it is
jointly continuous or uniformly delta-differentiable at the center.

With `l0=G_L(jL sigma)`, the approximation is

```
Qapp=G0[jH sigma-Bop l0],
gapp=ag*l0+bg*Qapp,
||g-gapp||infinity <=1000*r^4*||sigma||infinity.
```

At r=1/100 the absolute error is at most 10^-5 times the source supremum,
uniformly over the admitted delta and spatial momentum. The proof bounds
both the coupled feedback and the difference between G_H and G0. It
retains the actual delta-dependent observable map and both normalization
derivatives. This is a controlled nonlocal response approximation, not a
local-in-time derivative expansion or a general light-only EFT.

At fixed r this error bound does not shrink with delta. It is not a
relative error, a uniform derivative error, or an error budget capable of
resolving an O(delta) locked-cone effect. The bound does not establish
that the actual omitted response stays large either.

The source-free prepared sector has a different symplectically projected
source factor, `a^3/(2*Keff)`. Its center limit is 2/5; the complementary
part of the full diagonal is 8/5. This identifies a missing source
channel, but an impulse diagonal is not a low-frequency error lower
bound. Separately, nonzero H1_0(-r,r) pulses obey the Poincare RMS bound.
For the explicitly defined instantaneous stiffness proxy this gives
`Omega_sigma/m_proxy>6/13`. That only excludes a parametrically small
version of this particular compact-pulse hierarchy. The proxy is not a
rolling spectral gap or cutoff; persistent and approximately bandlimited
sources, longer preparation and correlated initial states remain
distinct research questions.

## What the global reciprocal calculation excludes

The [reciprocal theorem](../problems/P8/s6/matching/variable/global/reciprocal/FORMULATION.md)
is independent of the unpromoted G1 lapse ansatz. It uses the full
interaction polynomial `P=2(beta1+2*beta2*y+beta3*y^2)`, including clock
dependence of the coefficients, and the actual f metric equation.

Assume b=alpha/a, positive P and lapse c, no additional f null stress,
positive Einstein terms, a nondegenerate physical-g bounce and
`(a*a')'>=0` on the specified future domain. Then the bounce equation
forces c>y locally. An integrating factor preserves this inequality
globally, gives z=h/c>0 and z'>0 for u>0, and proves

```
integral from u0 to infinity of b*c du <=alpha/[z(u0)*a(u0)], u0>0.
```

Thus future f-null affine length and noncomoving f-timelike proper time
are finite in that domain, while physical-g causal geodesics are
complete. The quartic scale a=(1+u^2)^2 satisfies the required geometric
sign exactly. Time-reflected hypotheses give the corresponding past
statement. The theorem uses neither an assumed TT mass sign nor G1's
separate asymptotic mass calculation.

It does not prove a curvature singularity, common inextendibility,
incomplete comoving f observers, global scalar instability or a no-go
for nonreciprocal scales or extra f stress. Original P8 does not require
an added second metric to be complete merely because it appears in a
matching ansatz. Physical-g completeness remains intact, so this is not
an exclusion of the original C/D operator row.

## Remaining gates and verification

The controlled-source calculation now has an explicit causal memory
and omission budget. A local light-only matching claim would still need
an appropriate source/state class, actual mode hierarchy and a justified
local expansion with an error commensurate with the claimed effect.
The same parent must also support the required vacuum and canonical
matching construction; the original C/D operator dictionary, omitted
operators and loops, and the adopted finite-gravity dispersion remainder
cannot be inferred from these classical local results.

The reciprocal certificate pins 11 sources and the recursively rebuilt
S6.20 lineage. Its SHA-256 is
`c4feedf4a5137aafb02dacd97b57a0dc417ffaa6d2f9fb44d46cfbe11e2c7dc2`.
All 48 ordinary tests passed in 129.49 seconds. A separate fresh ordinary
CLI rebuild exited successfully and matched the frozen report exactly.
The 12-test independently authored audit reconstructs the full potential,
Einstein variations, integrating factor and geodesic bounds. There are
8 exact source/differential residuals and 16 invalid-input controls.

The frozen [forced-response certificate](../problems/P8/s6/matching/variable/response/forced/certificates/causal-physical-response.json)
has SHA-256
`0050320c15694cab2cc40aa8c3e82eaa345283a29d8347c0990e999c669ac755`.
It pins 18 sources, 26 exact residuals, 27 strict continuous margins,
162 independent nonnegative polynomial coefficients, 120 independent
Fraction calibration comparisons and 26 rejected-domain calls. All 72
ordinary tests passed in 147.85 seconds; a separate fresh ordinary
`--check` replay passed in 134.46 seconds. Both RNG seeds were zero,
with no watchdog timer or arithmetic adapter in these ordinary runs.
The 16-test independently authored root audit covers the literal action, moving-map jets,
both conserved TT polarizations, source normalization, Liouville and
Gauss identities, continuous pole estimate, independent resolvent error
chain and projection/RMS controls. The combined P8 regression passed
2,712 tests in 603.16 seconds, using the documented exact-arithmetic
adapter with 128 startup self-checks, 6,509 exact descents and 4,527
out-of-domain fallbacks. Both RNG seeds were zero and the faulthandler
plugin was disabled. Ruff, source hashes and whitespace checks passed.
The new working fixed-pulse phase and local-vacuum children were excluded
from this checkpoint's regression and commit manifest. No frozen ancestor
or unrelated P4/P9 work was changed by either certified child.
