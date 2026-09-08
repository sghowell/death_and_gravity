# Twice-vary the fixed counterterm before the dimensional limit

The common loop normalization is 1/(64 pi^2 epsilon_DR).
Use exactly the S6.53 scalar-coefficient continuation,

    m^4 C4(a_m,b_m)+m^2 R_D+2 a4V_4
      -m^2[(a_m-1)ga2+(b_m-1)gb2]
      -[(a_m-1)ga4+(b_m-1)gb4],

and add twice the S6.63 quadratic derivative mass pole at orders
two and four. C4 already contains the full potential and its
quadratic part; no separate potential quadratic pole is added.

For N=1+e n and a_hat=a exp(e zeta), the relative volume is

    N exp(D e zeta)
      =1+e(n+D zeta)+e^2(D n zeta+D^2 zeta^2/2)+O(e^3).

The physical expansion rate is (H+e zeta')/N. Generate its successive
proper-time derivatives with N^-1 partial_u. This construction
varies the curvature and the unit-normal coefficient tensors before
restricting to the clock. The tensors ga/gb are precisely the fixed
S6.53 tensors, including three times the spatial average, not a
replacement D-dimensional trace.

Quadratic mass deviations need only their first-order mass jets:
varying their geometric coefficient, volume, or a second mass jet
would be third order. In contrast, the terms linear in mass
deviations require the varied geometry, varied volume and both
second mass jets. All are retained.

Pure constant-coefficient box R terms are removed only as complete
compact-support Euler-null actions. The mass-weighted box R in
ga4/gb4 is retained. Curvature-squared invariants retain their full
D dependence; no four-dimensional Gauss-Bonnet identity is imposed
before variation and continuation.

## Pole and finite matching

Let Q_D be the quadratic density relative to the background a^D.
For output X its normalized Euler operator is

    E_X,D=sum_j[-(partial_u+D H)]^j partial Q_D/partial X^(j).

The separately computed radial Laurent poles agree with E_X,3 at
all three adiabatic orders, for both currents, on the actual clock.
No equality of arbitrary-D pole extensions is asserted.

Write D=3-2 epsilon_DR only after variation. The bare subtraction
-E_X,D/epsilon_DR contributes the finite term
+2 partial_D E_X,D at D=3. Add this to the radial finite part at
mu=m. This accounts for the volume's dimensional derivative as
well as the invariant and polarization derivatives.

The resulting matrix of normalized Euler-current operators is
formally self-adjoint in the background measure a^3 du. For an
entry sum_j c_j partial_u^j, its weighted adjoint has coefficients

    (K^dagger)_r
      =sum_(j>=r)(-1)^j binomial(j,r)
         (partial_u+3H)^(j-r)c_j.

All diagonal and mixed coefficients pass this identity. This
symmetry is a property of the local action Hessian; the full
retarded nonlocal response is not asserted to be time symmetric.

## Independent controls

At zero derivative order the result equals the Hessian of

    -N exp(3 zeta) m^4 B(a_m,b_m)

with the frozen finite potential B at mu=m. Thus the second mass
vertices and volume contacts are included without double counting.

Switching off all mass-deviation jets gives the independent
ordinary-Proca finite heat action

    (5/2)m^4+(5/3)m^2 R-4 a4sc_4.

Its two-current Hessian agrees at all three orders. Constant
spatial rescaling also gives K_NZ,0=-3 rho_local and
K_ZZ,0=9 p_local, with the same common normalization.

If the finite dimensional counterterm is omitted, the mixed
adjoint defect has bounce unit-mass fixtures -1600/81 at order
two and (-64/3,-512/27) in the zero/second derivative coefficients
at order four. These nonzero controls prevent a component-wise
radial finite assignment from masquerading as the matched action.
