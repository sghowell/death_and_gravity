# Complete original finite scalar Hessian, including the remainder

The source-pinned original four-dimensional finite density before64pi^2 is

5m^4/2 + 5m^2 R_old/3 - Weyl^2/30 + Euler/90 - R_old^2/18,

with m=1000. No finite coefficient is inferred merely from a cut or changed to simplify a factorization. Compact Euler and box-R variations vanish as actual variational boundaries; they are not pointwise zero densities and no infinite-history flux is silently discarded.

Let Ti=3wi-ci and Si,Wi be the first-curvature channels in FORMULATION.md. The mixed rescaled scalar curvature in the original exponential chart is

rDG=a^2 delta_D delta_G R_old
=24wD'wG'-8(wD'cG'+cD'wG')+4cD'cG'
 -20q wD wG+4q(cD wG+cG wD).

There is no freedom to replace the exponential second metric variation by a linear metric chart.

One independent derivation starts from the ADM extrinsic-curvature expression and the warped spatial three-metric scalar curvature. It gives the full nonlinear identity

a^2 R_old=6U+6w''+18h w'-2c''-6h c'
 +12w'^2-8w'c'+2c'^2
 +exp(-2(w-c))[-4w_xx-2w_x^2-4c_x w_x].

The test expands this separate formula through both Fourier legs and compares it to production's literal four-dimensional Christoffel contraction. It then reconstructs the complete local action, including the volume jet exp(3w-c).

The exact mixed R_old^2 density is

72SD SG+12U rDG+72U(TD SG+TG SD)+36U^2 TD TG.

The Einstein density is a^2[rDG+6(TD SG+TG SD)+6U TD TG], and the volume density is a^4 TD TG. Conformal invariance of the four-dimensional Weyl density, on a zero-Weyl background, gives its complete mixed contribution-(4/45)WD WG. Thus the curvature-square factor is-4SD SG-(4/45)WD WG, and the COMPLETE remaining local density is

Lrem=-(2/3)U rDG-4U(TD SG+TG SD)-2U^2 TD TG
 +(5m^2/3)a^2[rDG+6(TD SG+TG SD)+6U TD TG]
 +(5m^4/2)a^4 TD TG.

Every field jet of orderj is assigned weightj and k is assigned weight1. The exact polynomial contains total differential orders0,1,2 and no fourth-order remainder. Coefficient derivatives in a variational integration by parts do not increase that maximal field differential order. This is only an order statement: q terms need not define a bounded perturbation on the reference C_tH^r domain. No norm smallness follows from derivative counting.

As a nontrivial original-input check, setq=0 and wi=ci/3. The spatial amplitude is trace-free and its Frobenius norm squared is8ci^2/3. Dividing the mixed Hessian by8/3 gives

-cD''cG''/30 + (5m^2 a^2/6-U/3)cD'cG'.

This equals the actual S189 unit-Frobenius local tensor Hessian: the kinetic coefficient is2a^2 A, A=5m^2/12-R0/36, and R0=6U/a^2. Both coefficients and the full reconstruction are checked. This bridge confirms the finite local normalization; it does not identify the actual nonlocal curved response with the flat factors of Aref.
