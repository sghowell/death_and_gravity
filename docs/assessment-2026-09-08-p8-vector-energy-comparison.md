# P8: physical vector energy and a finite rolling-mode comparison

Original P8 remains open. The scoped P8(a) photon objective, original
linear CD/M1 classification and physical-frame V/G/B contract are
unchanged. This follows the [local curvature audit](assessment-2026-09-08-p8-local-vector-curvature.md).

## New result

[S6.50](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/state/FORMULATION.md)
derives the retained vector's actual contribution to the physical lapse
equation and isotropic pressure. It retains the temporal constraint,
canonical normalization and clock mass derivatives a_N=4/(9h),
b_N=28/(81h). These derivatives do not vanish when a=b=1 on the
clock. Varying before or after eliminating the temporal constraint
agrees; setting the coefficients to one before varying does not.

The zero-derivative energy pole independently agrees with C+C_N of
the local determinant, namely 3+20/(9h), not C=3. A direct spatial
dimensional Gamma integral fixes its sign and normalization. The
physical energy is not simply the canonical oscillator Hamiltonian.

For each actual transverse/longitudinal canonical rolling mode on
u in [-1/2,1/2], the checkpoint constructs a positive fourth-order
WKB reference W, with exact residual bounded by C/omega^4,
C=2000000, at m=m0*tau>=1000. Continuous rational coefficient
envelopes cover all momenta; no momentum grid or stationary metric
substitution is used. The reference has canonical Wronskian i.

Exact modes are prepared to match that reference and its derivative
at u=-1/2. Exact variation of constants preserves the CCR and bounds
the difference of the physical energy forms. The integral over all
R^3 comoving momenta is absolutely convergent and uniform in time.
Including two transverse and one longitudinal polarization gives

    abs(Delta rho)/(M^2/tau^2)<=12C/[(m0*tau)*(M*tau)^2],
    abs(Delta p)/(M^2/tau^2)<=108C/[5*(m0*tau)*(M*tau)^2].

At M*tau=10^12,m0*tau=1000 these are 2.4e-20 and 4.32e-20.
These are bounds on the exact-mode/reference differences, not on
either divergent unsubtracted integral or the total quantum correction.

## Remaining work and current continuation

The reference needs a complete subtraction and finite covariant
counterterm matching before the displayed number can be used in a
renormalized bounce equation. Finite-order Gaussian preparation is
not automatically an all-order Hadamard condition. Other field loops,
all-time state tails, corrected constraints/cones, interacting cutoff
and V/G/B remain separate work. No user intervention is needed now.

The next checkpoint derives a finite, explicitly prescribed adiabatic
subtraction and bounds its rational reference tail. Subsequent local
matching must retain dimension-dependent polarization/friction and
counterterm variations: a component-wise finite Laurent coefficient
need not be a covariantly renormalized stress component.

## Verification

The report pins 13 local sources and fully rebuilds S6.49 and its
ancestry. It checks 40 named exact identities comprising 45 scalar
entries, 17 continuous proof checks and 40 rejected inputs. Scientific
tests pass **15 tests in 75.95 seconds**. The ordinary suite passes
**45 tests in 334.93 seconds**, without the broad GCD adapter; the
separate seeded read-only CLI passes.

The full current P8 regression passes **4755 tests in 910.79 seconds**,
including every frozen checkpoint and the nine then-collected working
subtraction science tests. The independently controlled exact GCD
adapter passes 128 original tuple comparisons and records 5761 domain
fallbacks and 6509 exact descents under the per-test seed recipe.
No frozen checkpoint is excluded. The source manifest and missing,
extra or mutated report fields are checked. This is exact symbolic
verification plus written continuous proofs, not proof-assistant
formalization or external peer review. Original P8 is not closed.
