# S6.53: finite clock-mass matching in an explicit continuation

Status: the vector-only on-clock first metric variations are matched
and bounded in the named scalar-coefficient continuation. Original
P8 remains OPEN. This extends S6.52 without changing any frozen action,
physical frame, state preparation, finite potential or lapse jet.

## Fixed prescription and task

Retain the S6.42 vector action, S6.47 finite potential, S6.50 exact
Gaussian preparation, S6.51 full fourth-order subtraction and S6.52
dimensional radial calculation. Use the original timelike clock tube,
D=3-2epsilon_DR, scalar mass deviations delta_a=a_m-1 and
delta_b=b_m-1, and the unit clock normal. Both deviations vanish
along the clock; their lapse jets alpha,beta do not. Their first
homogeneous physical spatial variations vanish.

Specify the local counterterm density before metric variation. The
[prescription](notes/prescription.md) retains the full frozen scalar
coefficient C4=2b_m^2+a_m^(3/2)b_m^(1/2), independent of D, and the
ordinary covariant curvature counterterms. Its additional curvature
terms are linear in the mass deviations, with coefficients built
from Rg/12-5Ricci/6 and the fixed four-dimensional scalar curvature
pole stress plus g*boxR/6, evaluated in D+1 dimensions. Use three
times the spatial average, not an unannounced D-dimensional trace.
These choices explicitly fix the evanescent finite terms for this
first-variation problem. They are not unique matching data imposed
by an unknown UV parent.

## Established calculation and regulator limit

Vary before restricting to the clock and before taking the dimension
limit. The added mass terms have the required lapse variation and
zero first spatial variation. All three physical lapse poles match
S6.52. Subtracting the epsilon coefficient of the continued action
variation from the radial finite part preserves the frozen flat
finite term and reproduces the ordinary-Proca control when the lapse
mass jets are switched off. At mu=m the second-order local energy
coefficient is

    -(10+alpha-3beta)H^2-(2/3)beta*H'.

The fourth-order coefficient, the full dimensional scalar-tensor
calculation and actual-clock substitutions are retained exactly in
the source and certificate. Compactly supported local Euler variations
of boxR vanish; no infinite-time flux is erased.

The [regulator-limit proof](notes/dimension-limit.md) continues an
analytic pair of prepared exact modes, without conjugating complex D.
On |D-3|<=1/4 it proves a uniform integrable momentum envelope for
the exact-mode/reference difference and the reference subtraction
tail. The finite subtracted integral is holomorphic near D=3 and
tends to the physical S6.51 integral. This justifies combining that
integral with the dimensionally matched local coefficients.

The [continuous physical bounds](notes/bounds.md) apply on the full
interval u in [-1/2,1/2], all real comoving momenta, m0*tau>=1000
and the stated fixed subtraction prescription. At M*tau=10^12,
m0*tau=1000, each vector energy/pressure magnitude is below 10^-14
of M^2/tau^2. This is an absolute first-variation estimate with an
explicit error decomposition, not a sign or principal-symbol result.

## Boundaries

Only the retained Gaussian vector at fixed background and its first
homogeneous metric variations are covered. Higher mass-curvature
variations, all-order Hadamard admissibility, derivatives of the
nonlocal state term, the full clock equation, other field/higher
loops, unknown finite matching operators and the corrected bounce,
constraints, cones and interacting cutoff are not established.
The clock-mass energy exchanges with the clock equation and is not
separately conserved like ordinary Proca energy. Small unsigned
corrections do not protect a classically saturated matter cone.

The timelike projector prescription is not asserted to extend to
X=0. A common healthy Lorentz-invariant vacuum, finite-gravity Regge
remainder and the adopted V/G/B contract remain research requirements.
No user choice is required for the presently missing estimates.
