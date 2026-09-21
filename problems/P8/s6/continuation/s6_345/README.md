# S6.345: finite selected-box curvature coefficient

This checkpoint determines one known ordered-box loop coefficient in an
explicit covariant jet convention. It does not assign the independent
extra parent-theory curvature coefficient.

The entire off-shell flat degree6 polynomial is projected onto six fixed
Bose-symmetric covariant jet words. All220 flat monomials establish the
rank6 basis. The full96 original mass-ordered metric insertions, including
all connection terms, give a generic curvature difference with2430 exact
coefficient residuals:

    Delta M5_box,6 = chi_box*T/sqrt(kappa),
    chi_box = g^4*c_box(n)/(16*pi^2),
    c_box(n) = -4/315 integral_0^1
      t*(18+18t+92t^2+363t^3)/[1+(n-1)z]^5 dz,
    t=z(1-z).

The coefficient is strictly negative. Exact checks give
c_box(1)=-58/945 and lim n^2*c_box(n)=-2/105.
The unchanged original source satisfies
|chi_box|/A0<g^2*n/14400<10^187. This is a LARGE upper bound,
not perturbative smallness or a physical Taylor-truncation estimate.

The triangle in S344 uses a different explicitly stated local lift.
It must be converted, and every other sector treated consistently,
before an aggregate coefficient can be reported. Neither coefficient
is an extra term to add to the already complete S342 radiation.

The21 immutable source files comprise9 modules, two tests, this README,
FORMULATION and eight proof notes. The20-field native report contains
753 named exact checks,3206 scalar entries,43 gates,eight controls and
727 rejected inputs. External acceptance is recorded separately.
Original V/G/B/P8, independent parent matching, internal gravity,
inclusive probability, Regge, same-parent bounce and UV remain open;
scoped P8(a) is unchanged.

