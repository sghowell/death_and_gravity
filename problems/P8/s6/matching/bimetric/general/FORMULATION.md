# S6.4.HR: arbitrary-beta regular-flat background obstruction

For positive Einstein coefficients and arbitrary constant real HR
interaction parameters, a regular common spatially flat FLRW solution
cannot have a physical stationary slice with `H_g=0` and `D_g H_g>0`
under the matter assumptions below. The result includes the CD
nondegenerate bounce and a quantitative finite-window mismatch.

This extends the exact background screen, not the positive vacuum or
flat-source matching claims of the pinned beta1 parent. No vacuum health,
mass gap, weak coupling, or cosmological stability is assumed or concluded
for arbitrary parameters. General S6 matching and P8 remain open.

## Fixed parent and physical frame

Use the same P8 `+---` convention, `R_FLRW=-6(DH+2H^2)`, with

\[
S=-\frac G2\int\sqrt{|g|}R_g-\frac F2\int\sqrt{|f|}R_f
 -\mathfrak m^4\int\sqrt{|g|}\sum_{n=0}^4\beta_n e_n(\sqrt{g^{-1}f})
 +S_g[\Psi_g,g]+S_f[\Psi_f,f].
\]

`G=M_g^2>0`, `F=M_f^2>0`, `mathfrak m^4>0`; all five `beta_n` are arbitrary
constant real numbers. Zero interaction is included by zero coefficients.
The physical metric is the prescribed `g`, with no derivative metric
redefinition or additional old DHOST action.

Require a common flat FLRW foliation,

\[
g=N_g^2dt^2-a^2d\mathbf x^2,\qquad
f=N_f^2dt^2-b^2d\mathbf x^2,
\]

and smooth positive finite `a,b,N_g,N_f` on a neighborhood of the candidate
stationary slice. The square root has positive eigenvalues `(c,y,y,y)`,
`c=N_f/N_g`, `y=b/a`.

The primary case is minimally coupled, separately conserved classical NEC
matter only on `g`: `S_f=0` and `rho_g+p_g>=0`. A separate extension permits
an independent matter sector on `f`, satisfying `rho_f+p_f>=0` in that
metric and its own conservation equation. Positive-field-metric canonical
matter, with arbitrary potential, is an example in either sector.
The stress is homogeneous/isotropic as required by the FLRW background.
These are two separate matter actions, not one field coupled to both
metrics, and no cross-sector matter interaction is introduced.

## Full pointwise branch coverage

Define `P(y)=mathfrak m^4(beta1+2 beta2 y+beta3 y^2)`. Full lapse-first
variation and separate conservation give the undivided equation

\[
P(y)(N_g\dot b-N_f\dot a)=0.
\]

At `H_g=0` there are two exhaustive cases.

- If `P(y)\ne0` at the slice, continuity gives a neighborhood where it
  remains nonzero. The dynamical relation then gives
  `H_f=H_g/y` there and `D_fH_f=D_gH_g/(cy)` at the slice. The two null
  equations imply
  \[
  -2(G+Fy^2)D_gH_g=(\rho_g+p_g)+cy^3(\rho_f+p_f)\ge0.
  \]
- If `P(y)=0` at the slice, the `g` interaction null stress vanishes
  directly, giving `-2G D_gH_g=rho_g+p_g>=0`. No dynamical-branch relation,
  condition on `H_f`, or division by `P` is used in this case.

Simple, double and isolated roots, switching points, algebraic intervals,
and the identically zero polynomial are covered. Only pointwise
nondegenerate stationary bounces are excluded; this does not assert an
exclusion of every conceivable degenerate extremality.

## Evidence and limits

The [proof](notes/no-bounce.md) derives all four lapse/scale equations,
the root split, and inherited CD window bounds. A separate exact Fraction
Laurent-polynomial engine independently checks coefficients rather than
sampling parameters. The read-only [certificate](certificates/regular-flat-hr.json)
pins and replays S6.3.beta1 and its lineage; it does not alter them.
The [primary-source audit](notes/sources.md) distinguishes compatible known
identities from the specific deduction certified here, with no novelty claim.

This is not an exclusion of non-flat spatial curvature, non-common or
non-bidiagonal metric geometries, singular/other square-root branches,
quantum or NEC-violating matter, nonminimal/derivative matter, shared
doubly coupled fields, other physical-frame prescriptions, or general
bimetric modifications and UV completions. No arbitrary-beta vacuum
or curvature-squared matching coefficient is supplied by this theorem.
