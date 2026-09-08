# A fixed, selected-state scalar profile

This is a new action, not an alteration of a frozen counterterm or
an implicit change of quantum state. The adopted S6 contract allows
recorded on-tube corrections with bounds and imposes no naturalness
prior. This particular correction is deliberately state-tuned.

## Define the two coefficient functions once

Fix L=M*tau, m=m0*tau, the S6.55 all-order Cauchy data and the
S6.53 finite prescription on the original complete clock history.
For each k_com, use `p8_vector_hadamard.preparation.initial_data`
at u0=-1/2 and propagate the unchanged exact transverse and
longitudinal oscillator equations. Their normalized initial values
are v0=(2W_B)^(-1/2) and
v0'=(-W_B'/(2W_B)-iW_B)v0.

`state.integrands()` explicitly supplies the four physical quadratic
readouts, their fourth-order subtractions and their 2+1 polarization
sum. For either energy or pressure define the dimensionless profile
as L^-2 times the radial integral of that combined finite integrand,
plus the actual S6.53 local polynomial

    sum_{n=0}^2 C_n(u)m^(4-2n)/(64pi^2 L^2).

This definition retains the lapse-dependent mass readouts and all
finite matching terms. It is not ordinary Proca energy substituted
for the physical lapse observable, or normal ordering away its
finite local terms. The previous state and prescription are input
data, and their exact integral definition fixes the profiles uniquely.

## Smooth definition on the whole clock history

On every compact real-time interval containing u0, the background
coefficients, a^-1 and each finite-order derivative are bounded.
The exact oscillator has a unique smooth solution there. For any
fixed order N, sufficiently large k makes its local WKB reference
positive, uniformly on that compact interval. The S6.55 initial
data have the corresponding all-order asymptotics. Variation of
constants with the finite-order residual gives arbitrarily strong
high-frequency remainders by taking N large enough.

For a prescribed number of time derivatives, choose N still
larger to absorb the finite additional momentum powers in the
physical reconstruction and exact derivative rows. The local
adiabatic expansion cancels the fourth-order subtraction; the
remaining density has an integrable high-frequency tail.
At bounded k the massive oscillator data and readouts are smooth
and locally bounded; the regular k=0 limits are integrable with
the k^2 measure. Thus each differentiated finite integral is
continuous on the compact interval. Since the interval and
derivative order were arbitrary, both profiles are C-infinity
for every finite real u.

The same compact-interval large-k argument allows the frozen
dimensional continuation to be taken through the finite integral:
the mode initial data remain the same analytic continuation,
the bounded-k exact ODE is holomorphic in dimension, and the
subtracted high-k tail is uniformly integrable near dimension
three. No new finite part is introduced outside the original
window. This is a qualitative extension argument, not a uniform
all-time numerical estimate. The explicit C5 bounds used here
remain those on |u|<=1/2.

## Literal local action and first variations

Let rho_sigma,p_sigma denote those fixed functions, normalized by
M^2/tau^2. In the original physical matter frame add

    Delta P=(M^2/tau^2) w(x)
                [-p_sigma+(rho_sigma+p_sigma)(x+1)/2].

The smooth w(x) in notes/bounds.md is identically one on an open
neighborhood of the clock tube and zero near x=0. It gives an
explicit extension of this candidate, not a UV-completion claim.
Although the functions were designed using one state, they are
now fixed ordinary scalar coefficient functions of u=phi/tau.
Their evaluation need not be local in the *design calculation*;
the resulting action is local in its dynamical fields.

No vector or affine connection appears in Delta P. In particular,
the vector operator, its original-clock Cauchy data and its finite
prescription are unchanged. This makes the construction non-circular.
The addition is of the same formal loop order as the retained
Gaussian vector stress. It is not a universal rule for changing
renormalization when the state is changed.

For a lower scalar P(u,x), the clock energy and pressure are
rho_P=2x P_x-P and p_P=P at x=-1. Here they give exactly

    rho_Delta=-rho_sigma, p_Delta=-p_sigma.

At fixed metric, the clock Euler source is
P_u+2(P_x)'+6H P_x. It equals
rho_sigma'+3H(rho_sigma+p_sigma), the negative of S6.54's
actual vector Ward source. Therefore energy, pressure and clock
first variations cancel. The unchanged free-matter equation also
holds. The aligned vector mean source vanishes on the clock,
as before; the source-squared Gaussian term has zero first
variation there.

Consequently the original complete clock geometry is an exact
background solution of the reduced retained-vector Gaussian
semiclassical equations for this new action and its selected
state. Its completeness is geometric and unchanged. This is
not an assertion about the quantum equations of omitted sectors.

## Point-chart terms and the limit of the cancellation

Use the actual map a_phys=a_hat*exp(omega),
omega=-log[(h-1+N^-2)/h]/4, and a_hat=a*exp(v).
On the flat window, the added normalized density is
N exp(3omega+3v)[-p+(rho+p)(1-N^-2)/2].
Its linear coefficients are rho-3p/(2h) and -3p: these are
the negatives of the actual S6.57 lapse and scale forces.
Its quadratic coefficients, with rho,p held fixed, are

    n^2: Delta J=[(3/h-1)/2]rho
                    +[-1/2+9/(4h)-21/(8h^2)]p,
    nv: 3rho-9p/(2h),
    v^2: -9p/2.

The exact second differentiation is checked before imposing
the clock. No lapse/auxiliary velocity or spatial derivative
is added by this lower scalar. Its literal lapse-square
coefficient is bounded as in notes/bounds.md; this is not a
Dirac or principal-symbol analysis of a nonlocal quantum action.

Under subsequent metric or state variations, the profiles must
remain fixed. The actual vector covariance variation, contact
terms, renormalized retarded kernel and initial-state response
do not disappear. Another state leaves the difference between
its stress and the selected profiles. Full coupled cones,
light/heavy interactions, omitted loops, cutoff and V/G/B remain
open. In particular, a stable or causal quantum solution does
not follow merely from cancelling its background tadpoles.
