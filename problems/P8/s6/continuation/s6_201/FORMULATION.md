# Curved initial-prefactor removal with the complete finite correction

## Unchanged actual problem

Retain the same full parent, canonical mass1000 Proca sector,
kappa=10^800, all-order prepared CD state, compact unit slab and
fixed covariant prescription. D and Gamma are real smooth compact
spatial tracefree probes, with Gamma zero near the common preparation.
Their time supports may overlap.

At each original spherical computational regulator K define J_alpha
using the S6.197 constant-initial-alpha W8 readouts, and J_unit using
unit W8 in the same full memory and contact formula. Neither is asserted
to be a new exact state or a renormalized effective-action derivative.

## Complete finite correction

For each physical mode, |alpha0|^2=1+|beta0|^2. The source-pinned
initial mixing bound gives |beta0|<=B nu^-6, B=2e6 and nu>=1000.
Write b=|beta0|^2. Then b<1 and the exact memory weight increment is

    (1+b_k)(1+b_l)-1=b_k+b_l+b_k b_l
      <=2B^2(nu^-12+mu^-12).

The constant initial phases cancel. The full contact has increment b_k.
No momentum derivative or analyticity of the prepared state prefactor
is needed.

For M[f]^2=||f||L2^2+||grad_x f||L2^2, the full correction
R_alpha=J_alpha-J_unit has an absolute all-momentum limit with

    |R_alpha|<2e12 M[D]M[Gamma],
    |R_alpha-R_alpha,K|<1e17 M[D]M[Gamma]/K, K>=1000.

The full contact retains its one-mode regulator and all three modes;
the memory retains both created momenta and all nine pairs.
Both external-metric canonical factors give 8e-788 and 4e-783/K.

Strict displays apply for nonzero input norms. All corresponding
quantities vanish exactly if either test vanishes.

## Actual unit-reference decomposition

The finite-regulator identity

    J_actual-J_unit=(J_actual-J_alpha)+(J_alpha-J_unit)

and S6.197 give an actual-to-unit reference remainder below
3e23 M[D]M[Gamma], with error 2e27 M[D]M[Gamma]/K.
Unit-W8 amplitudes obey the same conservative time-analytic upper
bounds used in S6.198. Its exact six-step identity therefore gives

    J_actual,K=R_actual-unit,K+C_unit,K
               +sum_(j=0)^4 B_j,unit,K+F_unit,K.

Here F_unit is the last endpoint plus sixth bulk. The known finite
piece R_actual-unit+F_unit remains below 2e48 M[D]N61[Gamma],
with error 2e52 M[D]N61[Gamma]/K and canonical displays
8e-752 and 8e-748/K. N61 is the unchanged source time-H6 plus
spatial-H1 norm. The unknown full contact and first five equal-time
kernels are not included in these bounds.

## Boundary

Unit W8 has canonical Wronskian i but is not an exact mode-equation
solution or physical bisolution. The actual CD state is not replaced.
Joint complex spatial-momentum estimates, spatial subtraction and
original fixed covariant matching remain to be proved. No full mixed
inverse, nonlinear interacting background, stability, physical cutoff
or original V/G/B closure follows. Scoped P8(a) and A.20-A.23 are unchanged.
