# Complete physical ADM heavy finite action and fixed profile

Work with the unchanged S240 scalar prescription, not the different old vector finite coefficients. Put ell=log n and use R_old=-R_P8, so R0=6(H'+2H²) on the physical reference. Before the overall 1/(64 pi²), the entire four-dimensional finite heavy scalar action is

Cv + Ce R_old + Cr R_old² + Cw Weyl² + (ell/180) Euler,

with Cv=(3/2-ell)n², Ce=(ell-1)n/3, Cr=-ell/36 and Cw=-ell/60. Compact prepared variations make the Euler/divergence variation a boundary identity; no infinite-history flux is discarded. The S240 profile contains the exact constant -Cv/(64 pi² kappa0), which cancels the Cv action for every metric and clock field before either variation or constraint elimination.

## Full scalar geometry, including shift

Use physical proper-time ADM variables

ds²=N² dt²-a² exp(2zeta)[(dx+beta dt)²+dy²+dz²],
N=1+n_lapse, beta=partial_x B/a², u=t.

Introduce independent nilpotents eD,eG with opposite Fourier phases exp(-ikx),exp(ikx). Keep their full mixed coefficient, including metric, inverse metric, connection, Ricci tensor and scalar curvature second variations. No inverse Laplacian or reduced constraint is substituted.

The independent nonlinear ADM check uses rate=H+zeta_t-beta zeta_x, div=beta_x, K=(3rate-div)/N and KijKij=((rate-div)²+2rate²)/N². With
3R=-a^-2 exp(-2zeta)(4zeta_xx+2zeta_x²) and
Delta_gamma N=a^-2 exp(-2zeta)(N_xx+zeta_x N_x),

R_old=3R+KijKij+K²+(2/N)(partial_t-beta partial_x)K-2 Delta_gamma N/N.

Its four nilpotent coefficients agree with the literal four-dimensional Christoffel calculation. For either leg,

R1=6zeta''+24H zeta'-6H n_lapse'-12(H'+2H²)n_lapse
+2k² a^-2[n_lapse+2zeta+B'+2HB],

W=k²[n_lapse-zeta+B'-HB].

A direct Weyl component gives C0101=W/3. The scalar Weyl mixed density is 8 W_D W_G/(3a). The opposite-sign spatial phases have been retained. Under a complete infinitesimal time-coordinate change n_lapse=-T', zeta=-HT, B=T, the checks give R1=-T R0' and W=0. This is a geometry check, not a change in the fixed clock prescription.

Write V_m=delta_D delta_G sqrt(-g) and E_m=delta_D delta_G(sqrt(-g)R_old). The full Christoffel calculation gives E_m, including all its lower-order terms; it is not replaced by a principal symbol. Exactly,

V_m=a³[3(n_D zeta_G+n_G zeta_D)+9zeta_D zeta_G],

delta_D delta_G(sqrt(-g)R_old²)
=2R0 E_m-R0² V_m+2a³ R_D R_G.

The latter identity keeps the second curvature variation. E_m is serialized in full in the report.

## Fixed reference-profile terms and reference matching

Remove the already-cancelled constant from the full finite reference heat stress and write the remaining components before 64 pi² as

rho_bar=-2n(ell-1)H²+ell(6H²H'+2HH''-H'^2),

P_bar=2n(ell-1)(3H²+2H')/3
-ell(18H²H'+12HH''+9H'^2+2H''')/3.

They satisfy rho_bar'+3H(rho_bar+P_bar)=0 exactly. The corresponding profile is A_bar+B_bar(X-1), A_bar=-P_bar, B_bar=-(rho_bar+P_bar)/2. These reference functions are fixed under subsequent variations. Its ENTIRE mixed density is

a³[-(rho_bar+P_bar)n_D n_G
+3rho_bar(n_D zeta_G+n_G zeta_D)-9P_bar zeta_D zeta_G].

Consequently the complete matched finite local Hessian before 64 pi² is

(Ce+2Cr R0)E_m-Cr R0² V_m+2Cr a³ R_D R_G
+Cw 8W_D W_G/(3a)+the full displayed profile density.

The independent whole-action lapse and scale Euler variations vanish on the reference. Spatial linear terms are divergences; the homogeneous Weyl background is zero. Together with the exact Ward identity, this establishes reference matching of this local block, including the clock equation. It does not assert that the fixed state-profile piece is on shell by itself.

The remainder of the S240 fixed profile uses the exact state/subtraction reference integrals. Its full Hessian is retained and bounded separately in [bounds.md](bounds.md). The nonlocal determinant response is a different remaining contribution and is not discarded. The extra mass-one vacuum profile has zero clock quadratic jets.
