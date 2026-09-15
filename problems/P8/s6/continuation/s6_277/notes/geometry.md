# Physical volume, unchanged lapse and independent curvature reconstruction

The S262 source-pinned physical map is SPATIAL:

    g00=N^2 at zero shift,
    g_ij=-R^(-1/2) gamma_ij.

The physical lapse is N. Applying the spatial factor to g00 would be a
different metric and is not allowed. Both scalar and Maxwell ADM
coefficients independently select this same map and its volume factor.

For the diagonal Bianchi-I gamma of canonical.md,

    A1=e^alpha R^(-1/4)e^beta,
    A2=e^alpha R^(-1/4)e^(-beta),
    A3=e^alpha R^(-1/4),
    Vphys=e^(3alpha) R^(-3/4).

The normal congruence is comoving, unit timelike and geodesic because the
lapse is spatially homogeneous. With dproper=N du its expansion is

    theta=(1/N)d(log Vphys)/du
         =[3alpha_u-3(R_u+R_N N_u)/(4R)]/N.

All derivatives of R here belong to the FULL source. Only at the endpoint
use the licensed exact values R1,R_u0,R_N=-2. Along the actual full curve,
Vphys_tau/Vphys tends to-3J/2, while u_tau/tau tends to-DJ. Therefore

    tau theta tends to3/(2D).

The finite positive endpoint of Vphys does not imply finite expansion.

For independent geometric verification, construct all Christoffel symbols
of the general diagonal metric diag(N^2,-A1^2,-A2^2,-A3^2), with arbitrary
positive time functions. Its exact Ricci contraction in the original
(+,-,-,-) convention gives

    Rphys=-2 dtheta/dproper-theta^2-(H1^2+H2^2+H3^2),
    Hi=(1/N)d(log Ai)/du, theta=H1+H2+H3.

No FLRW formula is substituted for the anisotropic metric. The computation
retains the arbitrary homogeneous lapse and every directional scale.

The full canonical beta_u stays bounded at the endpoint. Hence
Hi=theta/3 plus a bounded anisotropic term, so their squared sum and
theta^2 grow only as O(tau^(-2)). In contrast,

    dtheta/dproper=3/(2D^2 J tau^3)+O(tau^(-2)),
    tau^3 Rphys tends to-3/(D^2 J).

The nonzero leading Ricci coefficient proves physical curvature divergence.
The local source and ODE have the derivative regularity needed to
differentiate these asymptotic quotients; one does not differentiate an
uncontrolled big-O bound. Equivalently, expand the smooth numerator and
denominator in tau first, divide by the simple zero u_tau, then
differentiate the resulting smooth-factor Laurent expression.

The lapse, directional scales and volume have finite positive limits,
but this geodesic reaches an unbounded curvature scalar after finite
proper time. A regular C2 physical-metric continuation through this
endpoint is impossible. No statement about weaker extensions, quantum
resolution or the validity of the literal EFT at arbitrarily large
curvature is inferred. In particular this is not a UV-completion no-go.
