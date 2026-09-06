# P8(a) A.12: absolute all-sampler QSEI on the local exact SEE solution

This checkpoint concerns the **new smooth exact solution of A.11**, not
the old off-shell A.8 metric. It derives a fresh bound from C1 potential
norms; it does not transfer A.9's higher-jet constants or its coefficient.

## Fixed field, state, geometry and prescription

Use one free real massless minimally coupled scalar, the original
positive Hadamard vacuum transported from its unchanged radiation past,
and the full A.11 local semiclassical solution. The ordinary classical
radiation remains `rho_rad=3 A^2/(kappa a_phys^4)`; there is no cosmological
term or added curvature-squared gravitational coupling. Retain the named
prescription `lambda=2 sqrt(2) A eta_star^2`, gamma=0 and physical
epsilon=1, delta=10^-14. In particular

    T0=A eta_star^2, kappa hbar=2880 pi^2 delta T0^2.

Here T0 is a physical normalization scale, not an asserted existence
lifetime. Set x=eta/eta_star-1, a_phys=A eta_star a, u=-a''/a,
h=a'/a and ds=a dx, where physical proper time is T0 times s up to
an irrelevant additive origin. Primes on a,u,h mean d/dx.

The old label y0=5/2 determines the matching slice x0. Let L=10^-10.
The actual unforced equation is established on

    I=(x0+L/2, x0+L).

The prescribed preparation source vanishes there. The new metric is
smooth and satisfies `2<=a<=3`, `1/3<=h<=1/2`, `|u|<=2*10^-12`.
On the entire original active history, whose duration is <=3,
`|u'|<=10^-5` and hence `|u|<=3*10^-5`. The past potential vanishes
smoothly at x=0. A.11 also proves
`||u'-ubar'||<72*10^-9` and common initial metric data.

No numerical u'' or u''' norm, higher curvature-jet cap, frequency gap,
or high-frequency vacuum replacement is assumed. The new reference
state is fixed by the original full history, not reset at x0.

## The absolute inequality

For every comoving timelike geodesic, every Hadamard target state omega
of the same scalar algebra on this smooth spacetime, and every real
proper-time sampler f supported compactly inside I,

    integral f(tau)^2 E_omega(tau) d_tau
       >= -2 hbar/(16 pi^2) integral |f''(tau)|^2 d_tau,
    E_omega=<T_UU>_omega-<T>_omega/2.

The same statement holds for real H2_0 samplers on compact subintervals
of I. The target need not be homogeneous, quasifree, zero-mean, close to
the reference, or itself a state solving SEE on this fixed metric.
The reference metric/state *does* solve SEE there. The finite scheme is
not varied when stating the absolute inequality.

The coefficient 2 is an explicit sufficient constant, **not a proved
optimum**. A rational root bound below 1+10^-8 is proved before rounding
its square upward to 2. The unusually small sampling interval matters.
There is no claim of a bound with this constant on the old full target
or on a longer exact-solution interval.

## Reference credit and scattering mechanism

The absolute reference credit is newly derived from the exact SEE and
fixed radiation, not imported from the old metric:

    E_reference >= hbar/(4199040 pi^2 A^4 eta_star^8) >0.

For normalized actual modes `w_k''+(k^2+u)w_k=0`, write
`w_k=e^(-ikx)b_k`, `q_k=u b_k`. The differentiated mode splits exactly
into forward and backward pieces. Integrating the latter once in its
past variable separates a local 1/k term, a real Fourier history of
u', and a 1/k remainder. Parseval in momentum controls the real history.
Together with full Parseval in sampling frequency this avoids demanding
pointwise high-frequency decay from unbounded higher potential jets.
The infrared modes are controlled directly by Volterra estimates.

## Focusing gate remains unfulfilled

The proper span available is at most `(3L/2) T0`. For the cosmological
comoving normals, the actual Ricci and extrinsic-curvature caps give a
strict lower bound on the FK index form which exceeds the available
initial contraction in either time direction. The QSEI-to-geometric
coefficient also has `Q2/tau_available^2 >=160000000` when the sufficient
coefficient 2 is used. These exact comparisons and the positive normal
Jacobian are proved in [notes/focusing.md](notes/focusing.md).

Thus this is a genuine quantitative QSEI application on a local exact
SEE solution, but not a cosmological incompleteness theorem. Other
hypersurfaces, boosted or null congruences, extended-time coverage,
realistic massive/interacting fields and semiclassical validity beyond
the prescribed mathematical model remain open. P8(a) and P8 are not closed.
