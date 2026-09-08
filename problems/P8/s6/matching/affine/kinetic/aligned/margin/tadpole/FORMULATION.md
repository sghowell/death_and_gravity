# P8-S6.60: selected-state clock-preserving scalar-profile candidate

Date: 2026-09-08. Original P8 remains OPEN.

This is a new named action under the adopted S6 permission for
on-tube corrections with an explicit error budget. It does not
overwrite S6.56 or its fixed-source response, change the S6.55
vector state or revise the S6.53 finite matching prescription.

## Literal candidate

Let rho_sigma(u),p_sigma(u) be the actual, dimensionless finite
stress profiles of that fixed predecessor state on the complete
original clock history, normalized by M^2/tau^2. They are fixed
c-number functions defined by the exact Cauchy problem and
convergent subtracted mode integrals in `state.py`.

To S6.56 add the Lorentzian lower-scalar density

    Delta P=(M^2/tau^2) w(x)
                [-p_sigma(u)+(rho_sigma(u)+p_sigma(u))(x+1)/2],
    u=phi/tau, x=g^{mu nu} partial_mu phi partial_nu phi,
    w(x)=1-chi(64(x+1)^2),

where chi is exactly the S6.55 smooth step, zero for s<=1 and
one for s>=2. Thus w=1 throughout |x+1|<=1/10, and w=0 in
an open neighborhood of the vacuum x=0. All other terms are
unchanged. The addition is assigned the same formal loop order
as the selected Gaussian vector stress.

This is deliberately finely tuned to one selected state and history.
The profiles are fixed before varying the new action; they are
not recomputed on a different metric or in a different state.
It is a physical scalar-action change, not a state-dependent
renormalization prescription or a universal cancellation of loops.

## Acceptance gates

1. Derive the retained vector's off-clock principal speed and the
   sign failure on the particular S6.57 approximate initial metric.
2. Give the exact, non-circular profile definition from the fixed
   vector Cauchy data, physical readouts and finite matching terms.
3. Prove local smoothness over the whole real clock history and
   quantitative mixed derivative bounds through total order five
   on |u|<=1/2, including the smooth clock-norm transition.
4. Verify cancellation of the selected energy, pressure and Ward
   clock source and derive the actual point-chart quadratic terms.
5. Record the selected-sector and perturbation boundary: first
   variations cancel, but the causal second variation is not removed.
6. Pin all sources, fully rebuild S6.59, reject unsupported inputs
   and pass the ordinary and full P8 regression suites.

## Scope

The original complete clock geometry and free-matter history solve
the reduced retained-vector Gaussian semiclassical background
equations for this new action and its selected state. Uniform small
coefficient bounds are claimed only on the declared bounce window;
global smooth definition is not a uniform all-time error estimate.

At M*tau=10^24,m0*tau=1000, all mixed (u,x) derivatives of the
dimensionless addition through total order five are below 10^-18
on |u|<=1/2 and every real x. Its clock-tube P_xx is zero, and
the literal point-chart lapse-square change is bounded by 15eta0/8.
No local field coefficient changes inside the open zero-gradient region.
If a predecessor vacuum lies there, its local jets are unchanged; this
does not establish vacuum existence or health.

No stability, full coupled quantum cone, interacting cutoff, complete
quantum solution, common UV parent, finite-gravity Regge or V/G/B
verdict is implied. Other quantized sectors and higher loops remain
outside the selected Gaussian calculation. Original P8 is not closed.
