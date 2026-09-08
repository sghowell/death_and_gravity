# S6.48: vector mass-insertion pole and finite subtracted loop bound

Certificate-gated, limited quantum component. Original P8 OPEN.
Keep S6.42 and the S6.47 Euclidean/dimensional conventions unchanged.
Calculate the vector-only two-insertion contribution from a varying
mass tensor at the isotropic on-clock point, with the physical metric
held fixed. Retain the full Proca propagator, including its longitudinal
numerator. Derive the Feynman-parameter and tensor-integral pole at
general external momentum, not just its zero-momentum potential.

Use the actual first coefficient variation A=diag(a_N,b_N,b_N,b_N),
a_N=4/(9h), b_N=28/(81h). The external r parameter labels this mass
variation. It is not a free physical lapse after all gravitational
constraints, nor a claim to the entire quantum scalar two-point function.
Separate the two-insertion contribution from the second mass variation's
local tadpole, metric/complement/light/gravity contributions and finite
matching counterterms. A divergent derivative coefficient is not a
finite remainder, a physical extra pole/cutoff or a UV no-go theorem.

Also derive the finite two-insertion loop remainder after Taylor
subtraction through fourth external-momentum order. Keep its explicit
logarithmic form factor, complex momentum domain and local-subtraction
boundary. Bound this particular loop contribution below the massive
threshold; do not promote it to a full curved/in-in or EFT remainder.

## Declared result and domain

The quadratic functional convention is the coefficient of r(-k)r(k),
with no additional implicit factor of one half. The pole is
-[F0+F2+F4]/(64*pi^2*epsilon_DR), with epsilon_DR=(4-d)/2 and
Fj homogeneous of external-momentum degree j. The full numerator,
angular moments, radial Laurent coefficients, parameter integral,
actual mass direction and independent isotropic reduction must agree.
On |u|<=1/2, h=(1+u^2)^3 and nonzero real Euclidean k,
(k^2)^2/125 < F4 < (k^2)^2/25. This is a derivative counterterm
coefficient, not a finite error or a physical propagating mode.

For the finite loop, subtract its momentum Taylor polynomial through
degree four. Let kappa=sum_mu |k_mu|^2 be the Hermitian norm of complex
Euclidean momentum, m0^2=M^2/zeta_physical, and rho=kappa/m0^2<=1.
The analytic logarithmic remainder obeys

    |R_loop(k)| <= m0^4/(64*pi^2*h^2)
                    * (4852/229635)*rho^3/(1-rho/4)
                <= kappa^3/(1920*pi^2*h^2*m0^2).

The last bound is strict for nonzero kappa; the zero-momentum remainder
is zero. With Lp=M*tau, R=m0*tau and Q=kappa*tau^2 in [0,R^2],
the ratio of this quadratic kernel to M^2/tau^2 is at most
Q^3/(17280*Lp^2*R^2). This uses a reference density, not a positive
lower bound on a bounce observable. A band-limited source-functional
bound additionally carries its squared Fourier/L2 norm.

Local finite counterterms through fourth momentum order disappear
under the specified subtraction. Unknown higher local matching
operators do not. No finite counterterm or original model is changed.
The complex ball permits low-momentum stationary analytic continuation;
it does not transfer the result to an evolving curved or in-in state.

## Verification gate

Pin and fully rebuild S6.47 and its ancestry. Pin every local source,
proof, test and formulation file. Require exact residuals, a written
continuous complex-log bound, strict scale/domain rejection, missing
longitudinal and incomplete second-variation countercontrols, and
read-only report replay with mutation rejection. This is not
proof-assistant formalization or a V/G/B admissibility certificate.
