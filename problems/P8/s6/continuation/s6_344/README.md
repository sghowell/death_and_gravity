# S6.344: finite selected-triangle curvature coefficient

This checkpoint computes a known finite matter-loop coefficient, not the
independent added parent-theory curvature matching of S336.

Relative to the explicitly specified scalar labeled-box lift, the full
selected triangle first differs at six derivatives. All three literal
internal line insertions give the generic simplex identity

    Delta_TT,6 = -8*x^2*y^2*z^2*Q / M^4,
    M=x+y+n*z,  x+y+z=1,
    Q=epsilon((k.R1)R2-(k.R2)R1, (k.R1)R2-(k.R2)R1).

The lower degree2 and degree4 differences vanish identically. Including
the original outer heavy branch, all24 scalar labels and the loop factor,

    chi_triangle = g^2*(C+g^2/n)*c_triangle(n)/(16*pi^2),
    c_triangle(n) = -4/15 integral_0^1 z^2*(1-z)^5/[1+(n-1)z]^4 dz.

The primitive moment is strictly negative, c_triangle(1)=-1/630, and
n^3*c_triangle(n)->-4/45. For the unchanged original source the dressed
coefficient is positive and |chi_triangle|/A0<10^-206.

The result is an off-shell analytic-origin derivative coefficient in an
explicit comparison convention. It is not a physical above-threshold
Taylor approximation, a truncation-error bound, or a new term to add to
the already complete S342 radiation. The independent extra parent chi,
boxes and aggregate matching, internal gravity, inclusive probability,
original V/G/B/P8 and scoped P8(a) retain their previous statuses.

The21 immutable source files include9 Python modules, two test files,
this README, FORMULATION and eight proof notes. The native certificate
has20 fields,514 named exact checks,536 scalar entries,42 proof gates,
eight controls and670 rejected inputs. External ordinary/CLI/full-suite
acceptance is recorded separately, never inside this source manifest.

