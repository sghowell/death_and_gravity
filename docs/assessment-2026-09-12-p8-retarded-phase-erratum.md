# P8 continuation: retarded pair-phase correction

Original P8 remains OPEN. Certification of the **physical finite-matching identification in S6.212 and the assembled-current identification in S6.213 is withdrawn pending a corrected successor**. Their frozen artifacts remain unchanged and reproducible; passing their tests does not resolve this scientific error.

## Independent finding

[The standalone exact audit](../scripts/p8_retarded_phase_audit.py) constructs quadratic observables directly from finite Fock matrices and independently from the real symplectic covariance tangent. Nine source/readout combinations agree with one another and with the physical pair orientation. The fixture contains mixed position/momentum vertices, so it tests complex stripped amplitudes rather than only the real flat case.

For annihilation pair coefficients b_D and b_G (with quadratic Hamiltonians normalized as Z^T M Z/2), the vacuum identity is

    i <[H_D(t),H_G(s)]> = +Im(conjugate(b_D(t)) b_G(s))
                        = -Im(b_D(t) conjugate(b_G(s))).

Writing b_X = a_X exp(-i Theta), the detector-sharp/source-annihilation form is therefore

    +Im[conjugate(a_D(t)) exp(+i Theta(t))
         integral_(t0)^t exp(-i Theta(s)) a_G(s) Gamma(s) ds].

The S195 canonical tangent and Kubo sign are consistent with this identity. The S197 independent Fock comparison also explicitly conjugates mode products to obtain creation amplitudes.

The S198 negative-phase integration identity is correct **for creation-pair amplitudes**, as its formulation specifies. The error is its later use with the opposite amplitude branch. In S207's curved extraction and S212's `jets.frequency`, the source momentum is proportional to **-i W minus the real amplitude derivative**, hence is an annihilation amplitude. The detector is its Schwarz partner, but the extraction retains the creation-amplitude coefficient and minus-imaginary-part rule.

For arbitrary stripped complex amplitudes, the incorrectly combined convention differs from the physical current by a nonzero cosine term. The standalone audit includes an exact nonzero counterexample; this is not a choice of Fourier sign or a new finite counterterm.

## Correct retarded identity

For the detector-sharp/source-annihilation convention, define

    Kplus(t)=exp(+i Theta(t)) integral_(t0)^t exp(-i Theta(s)) b(s) ds,
    Omega=Theta', g=1/Omega, Lb=(g b)'.

Then

    Kplus' - i Omega Kplus = b,
    Kplus = sum_(j=0)^(n-1) i(-i)^j g L^j b
            +(-i)^n exp(+i Theta(t)) integral exp(-i Theta(s)) L^n b(s) ds.

The current is **+Im** after contraction with the detector-sharp amplitude. The same initial germ removes the same lower endpoints. All six ODE identities are independently checked without importing a frozen endpoint helper.

Consequently, with the frozen stripped amplitude product held fixed,

    B_j(correct physical current) = (-1)^j B_j(frozen extraction).

Even endpoints are unchanged. Odd endpoints reverse sign. The bulk phase and imaginary-part convention must be corrected consistently too; adding a local counterterm to force a match is not a valid repair.

## What exposed the issue

Private full-ADM scalar work retained the trace-dependent longitudinal constraint and the complete metric volume. Literal finite vertices, integer-dimensional field-strength geometry and independent full four-order WKB expansions agreed with the algebraic scalar tables. Yet those tables did not match the original covariant curvature pole.

For unnormalized trace I and longitudinal tracefree scalar S=diag(-1,-1,2), let Delta be the original pole Hessian minus sixteen times the frozen-orientation normalized logarithmic density. At a=(1+t^2)^2,

    Delta_tt = 2 a'' p^2 Gamma,
    Delta_tS = 2 a' p^2 Gamma',
    Delta_St = -2 p^2 (a'' Gamma + a' Gamma').

The sign reversal of the calculated odd endpoints reproduces all three discrepancies exactly. This diagnosis preceded applying the correction. A fresh corrected calculation then matched all three independently constructed full-volume pole Hessians exactly (60.70 seconds). No scalar current is claimed complete by this calculation.

The inherited tracefree dimensional coefficient gives an explicit nonzero finite correction. In the fixed invariant basis T=tr(D Gamma), V=e^T D Gamma e, W=(e^T D e)(e^T Gamma e), the only nonzero odd spatial coefficient is j=1, source jet 0, radial degree 3. Its physical value is zero; its dimension derivative is

    p^2 (3t^2+1) (-31/420, 1/14, 0).

Therefore the corrected minus frozen S212 local finite term is

    delta Ffinite = p^2 (3t^2+1)/(420 pi^2) (-31 T+30 V) Gamma_0.

Here T,V denote contractions of the fixed detector and source tensor directions; Gamma_0 is the independent scalar time profile used in the coefficient convention. There is no correction to its first- or second-time-jet coefficient. The physical pole, lower-band power term and prescribed curvature counterterm are unchanged. This explicit correction alone does not repair the remaining finite-endpoint/bulk decomposition of S213.

## Scope and disposition

- S195's canonical tangent/Kubo identity, S197's actual-state comparison and S198's abstract creation-amplitude integration identity are not refuted by this finding.
- S207's identification of its complex annihilation-branch extraction with that creation-branch template needs correction. Every descendant using those oriented coefficients must use the corrected branch.
- The physical tracefree logarithmic odd spatial difference vanishes, but its dimension derivative does not. Thus the successful tracefree physical pole check did **not** validate S212's finite part. S212's physical finite-part identification is withdrawn.
- S213's old assembled-current and anchored-regulator identification cannot be used as a completed physical response pending consistent correction of finite endpoints, bulk, finite matching and cutoff subtraction. Its numerical bounds are not a substitute for that identification.
- Absolute-value analytic bounds are insensitive to a unit phase/sign, but their application to the corrected decomposition must be stated and checked, not silently assumed.
- S214's new local constrained ADM vertices and metric-chart contacts do not use the erroneous endpoint extraction. Their independent tests and 42063-test regression passed. Their inherited status table is historical, not a repair of S212/S213.
- S215's prepared Ward identities are conditional geometric identities. Its claim to a known tracefree input depends on the withdrawn S213 identification; that dependency is not currently certified.
- Scoped P8(a), the A.20-A.23 gates, and the original open V/G/B obligations are unchanged. This is not a new user-intervention blocker.

Frozen source, tests and certificate bytes are preserved as historical evidence. This notice supersedes contrary current-status language in the S212/S213 publication audits; it does not rewrite those historical validation records. Work continues on an explicitly corrected successor and stronger independent regression tests.
