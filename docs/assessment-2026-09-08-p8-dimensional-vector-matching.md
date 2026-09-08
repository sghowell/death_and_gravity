# P8: dimensional vector poles and finite local matching controls

Original P8 remains open. The scoped P8(a) photon objective, original
linear classification and physical-frame V/G/B contract are unchanged.
This follows the [finite subtraction audit](assessment-2026-09-08-p8-finite-vector-subtraction.md).

## New result

[S6.52](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/dimensional/FORMULATION.md)
continues the retained-vector action and mode calculation to D spatial
dimensions. Direct action variation retains the D-dependent canonical
rates and D-1 transverse polarizations. Twenty reductions recover the
actual frozen four-dimensional WKB and subtraction functions.

The radial Gamma calculation retains the dimension dependence through
the finite term. Taking the polarization count to four dimensions too
early loses 2+4b_N in the finite zeroth-order lapse coefficient. Its
correct value preserves the already frozen S6.47 finite potential and
N derivative. The actual order-zero, two and four lapse poles are
independently reconstructed from covariant one-form/scalar-gradient
coincidence coefficients; the constrained scalar-gradient term matters.

Writing alpha=a_N,beta=b_N and E=H'^2-6H^2H'-2HH'', the extra
clock-mass derivative poles, in the stated common normalization, are

    extra_P2=-(3H^2+4H')alpha-(9H^2+2H')beta,
    extra_P4=(beta-alpha)E/2+(alpha+beta)boxR/12.

These agree with the independent radial calculation. Only local
dimensional poles are used in the scalar EOM/trace argument; finite
anomalies are not set to zero.

For the separate ordinary-Proca control alpha=beta=0, the raw radial
finite parts fail the four-dimensional conservation diagnostic. Varying
the dimensional covariant counterterm action before taking the limit
supplies the finite terms that restore it. The resulting local energy
coefficients, with respective m^4,m^2,1 prefactors and common 1/(64pi^2),
are

    3ell-5/2,
    H^2(6ell-10),
    (ell+2)(6H^2H'+2HH''-H'^2), ell=log(m^2/mu^2).

An independent finite local heat-action variation agrees. The ordinary
control is not substituted for the actual clock-dependent finite energy.
Compactly supported local variations of boxR vanish; no noncompact
infinite-time flux is assumed to vanish.

## Remaining work

The actual finite clock-mass curvature counterterms still need an
explicit continuation preserving the frozen potential and lapse jets.
Different evanescent continuations can agree on a four-dimensional
pole but shift finite terms. The next candidate prescription is
recorded separately, with continuous local bounds under development.
Before combining its local terms with the finite mode result, a
uniform bound must justify the regulator-dimension limit of the
subtracted mode integral.

State admissibility, full scalar/quantum variations, all-time control,
other field/higher loops, corrected constraints/cones, interacting
cutoff and V/G/B remain open. No old action or finite potential has
changed. There is no user-intervention blocker at present; the missing
matching estimates are ongoing research.

## Verification

The report pins 14 local sources and fully rebuilds S6.51 and its
ancestry. It checks 71 named exact scalar identities and 33 rejected
inputs. Scientific tests pass **10 tests in 77.90 seconds**. The
ordinary suite passes **41 tests in 407.09 seconds**, without the
broad GCD adapter; the separate seeded read-only CLI passes.

The full P8 regression through S6.52 passes **4829 tests in 898.49
seconds**, without excluding any frozen checkpoint. The exact GCD
adapter passes 128 original tuple comparisons and records 5763 domain
fallbacks and 6509 exact descents under the per-test seed recipe.
Source manifests and missing, extra or mutated report fields are
checked. This is exact symbolic verification with written local proofs,
not proof-assistant formalization or external peer review. Original
P8 is not finished or closed.
