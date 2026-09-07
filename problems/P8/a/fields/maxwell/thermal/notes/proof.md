# Actual physical thermal state, exact branch and continuous C3 estimate

## 1. Construct the physical state, not a free radiation constant

Take the usual two transverse photon oscillators in flat conformal
coordinates (eta,X), with k=|k_vector|>0 and transverse polarization
projector sum_lambda e_i^lambda conjugate(e_j^lambda)
=delta_ij-k_i k_j/k^2. Their field amplitude is sqrt(hbar/(2k)) and
their occupation covariance is

    <a_dagger(k,lambda)a(k',lambda')>
       =(2pi)^3 delta^3(k-k') delta_lambda,lambda' n(k),
    n(k)=1/(exp(b_T k)-1),
    <a(k,lambda)a_dagger(k',lambda')>
       =(2pi)^3 delta^3(k-k') delta_lambda,lambda' [1+n(k)],
    <aa>=0.

For every positive b_T these nonnegative oscillator covariances define
a positive quasifree state on the physical field-strength algebra and
retain the same vacuum commutator. Gauge-potential representations are
unnecessary for this construction. The smooth thermal-minus-vacuum
two-point function of the physical field strength can be checked
directly: at large k, n(k) and all polynomially weighted integrands
decay exponentially. At small k, n(k)=O(k^-1); even the potential
amplitude k^-1 with radial measure k^2 is locally integrable, and the
field-strength derivatives improve the power. Further time or spatial
derivatives only add nonnegative powers of k. The bounded transverse
projector causes no angular divergence. Differentiation under the
integral is therefore justified to every order on compact sets.

The field-strength vacuum two-point function is Hadamard, and adding
this smooth thermal difference preserves that microlocal property.
The four-dimensional conformal Maxwell algebra map of A.16 transports
the state to any smooth positive flat-FLRW strip, preserving positivity,
the commutator and the Hadamard wavefront set. The physical E and B
components acquire a^-2 factors, and the thermal stress difference
therefore has physical a^-4 weight. This is an actual fixed conformal
thermal state on the resulting metric, not a sequence of time-dependent
instantaneous vacuum choices or a covariance frozen on another metric.

The actual Minkowski thermal density above vacuum is

    Q=2 hbar integral d^3k/(2pi)^3 k n(k)
      =hbar/pi^2 integral_0^infinity k^3/(exp(b_T k)-1) dk
      =hbar*pi^2/(15 b_T^4).                               (1)

For the last identity expand the nonnegative Bose factor into its
exponential series. Tonelli's theorem and
integral k^3 exp(-n b_T k) dk=6/(n b_T)^4 give
6 zeta(4)=pi^4/15. Isotropy gives pressure Q/3. The coefficient counts
two physical photon polarizations, not one scalar oscillator.

The length b_T is a conformal-time inverse frequency. In c=1 units the
locally measured physical inverse frequency is a b_T, so
k_B T_physical=hbar/(a b_T). Substitution into the ordinary physical
blackbody formula pi^2(k_B T)^4/(15 hbar^3) recovers (1)/a^4. Calling
b_T simply 1/(k_B T) while retaining explicit hbar is dimensionally
incorrect. Under a conformal-coordinate rescaling a->c a and
b_T->b_T/c, Q->c^4 Q and Q/a^4 is unchanged. Covariant T_etaeta/a^2,
not T_etaeta itself, is the physical energy density.

## 2. Prescribed reference anomaly and exact SEE components

The photon finite prescription is explicitly beta_M=0; Lambda=0 and
all additional sources or independent nonzero gravitational
curvature-squared terms are absent in this named example. These are
new explicit specializations inside A.16's already stated family, not
changes to A.18 or an identification with scalar gamma. A finite Weyl
variation vanishes in flat FLRW, but a nonzero finite R^2 term would
change the present dynamics and cannot be omitted without this choice.

Set r=31 hbar/(480pi^2). A.16's actual conformal vacuum and the
traceless thermal difference give, with physical positive density,

    rho=Q/a^4+rH^4,
    pressure=Q/(3a^4)-rH^4-(4r/3)H^2 Hdot,
    EED=(rho+3 pressure)/2=Q/a^4-r(H^4+2H^2 Hdot),
    trace_FK=rho-3 pressure=4r(H^4+H^2 Hdot).                (2)

Direct differentiation proves rho_dot+3H(rho+pressure)=0. Neither the
thermal charge nor the vacuum term is chosen from conservation alone:
the occupation state fixes Q, and A.16's conformal reference fixes rH^4.
The trace is not zero, despite the vanishing trace of the stress
difference. Omitting rH^4 would replace the exact equation by the
classical radiation approximation and change the initial constraint.

The FK metric is (+---), R=6(Hdot+2H^2), G_00=-3H^2 and the actual SEE
is G_FK=-kappa T. Its two physical diagonal components are

    3H^2=kappa rho,
    2Hdot+3H^2=-kappa pressure.                             (3)

Let b=kappa r/3. The first equation is

    H^2-bH^4=kappa Q/(3a^4).                               (4)

After substituting (4), the second equation is exactly

    2(1-2bH^2)Hdot+4H^2(1-bH^2)=0.                        (5)

Conservation and (4) imply (5) when H!=0, and the module independently
checks the direct substitution so no pressure equation is lost.
This is the full exact equation in the named prescription, not a
truncation in bH^2 or an order-reduced surrogate for another RSET.

## 3. Parameter and normal-clock match

Choose tau>0 with 0<delta=kappa hbar/(8pi^2 tau^2)<=10^-8 and set
lambda=b/tau^2=31delta/180. Use s=T_0-T for the contracting normal in
an expansion-oriented model; x=s/tau and y=-tau H_normal=+tau H_expansion.
Then (5) becomes

    y'=F(y)=2y^2(1-lambda y^2)/(1-2lambda y^2).              (6)

The scale factor obeys d_x log a=-y. Fix a(0)=1 and y(0)=2. Equation
(4) fixes Q=12(1-4lambda)/(kappa tau^2), which is positive on the named
parameter range. Choosing b_T=[hbar*pi^2/(15Q)]^(1/4) in the state of
section 1 realizes this precise Q. No arbitrary homogeneous radiation
constant or separate classical fluid is assumed. The relation between
delta, kappa, hbar and tau is retained; they are not independently
varied in a way that violates the source match.

For fixed positive delta, (6) is analytic and locally Lipschitz in y
on the low interval 0<y<ycrit=1/sqrt(2lambda). On this interval its
derivative is positive. The first integral of the metric is

    a^4=4(1-4lambda)/[y^2(1-lambda y^2)]>0.                (7)

Differentiating (7) with (6) gives d_x log(a^4)=-4y and its value at
y=2 is one, so it supplies the actual positive scale factor on the
same solution. The implicit clock is

    Phi(y)=1/(2y)+(sqrt(lambda)/2)atanh(sqrt(lambda)y),
    Phi'(y)F(y)=-1,
    x=Phi(2)-Phi(y).                                      (8)

It is monotone and single-valued on the chosen low branch. The use of
principal real atanh is valid because sqrt(lambda)y<1/sqrt(2)<1.
The high branch, even where the algebraic a^4 remains positive, is
not selected by these initial data or by the monotone low-branch clock.
The ODE does not permit a smooth switch through its critical denominator.

## 4. Actual short interval and the zeroth comparison

Put r0=1/100, not to be confused with the anomaly coefficient r in (2).
The named lambda is less than 1/16. On the provisional backward box
1<=y<=2, the denominator in (6) is at least 1/2, so

    0<F(y)<=4y^2.

Monotonicity ensures y<=2 to the past. Let v=1/y. Then v'>=-4, so for
x in [-r0,0], v(x)<=1/2+4r0=27/50. Hence y(x)>=50/27>1. The solution
cannot leave this compact box on that interval; standard ODE
continuation proves existence and smoothness on a neighborhood of its
closed endpoints. No small-residual-to-solution inference is involved.

For the anchored radiation reference y_rad=2/(1-4x), its reciprocal is
v_rad=1/2-2x. The exact difference equation is

    v'=-2-2lambda y^2/(1-2lambda y^2),
    v(x)-v_rad(x)=integral_x^0 2lambda y(t)^2/
                                      (1-2lambda y(t)^2) dt.

Therefore v>=v_rad>=1/2 and

    0<=y_rad-y <=4(v-v_rad)
       <=32r0 lambda/(1-8lambda)
       <=(16/25)lambda.                                  (9)

This includes every x of the interval and the anchor y(0)=y_rad(0)=2.
It is not a first-order expansion with an unknown remainder.

## 5. Three derivative comparisons with exact response polynomials

Write z=lambda y^2 and D=1-2z. On the backward box 0<=z<=1/4. Along
the exact solution y'=F, y''=F F_y, y'''=F(F F_y)_y. Their deviations
from the corresponding radiation polynomials evaluated at the same y
are exactly

    F-2y^2 =2lambda y^4/D,
    F F_y-8y^3 =8lambda y^5(6z^2-8z+3)/D^3,
    F(F F_y)_y-48y^4
      =16lambda y^6(6z^2-5z+2)(14z^2-22z+9)/D^5.          (10)

On [0,1/4] all three quadratic factors are decreasing and positive;
their lower endpoint-in-value numbers at z=1/4 are 11/8,9/8,35/8,
and their upper values at zero are 3,2,9. Use y<=2 and D>=1/2 to
bound the three absolute corrections by

    (64,6144,589824)lambda.

On 0<=y,y_rad<=2 the radiation polynomials 2y^2,8y^3,48y^4 have
Lipschitz constants (8,96,1536). Adding these constants times (9)
gives the complete continuous C3 error vector

    ||d_x^j(y-y_rad)||<=lambda C_j,
    C=(16/25,1728/25,155136/25,14770176/25), j=0,...,3.     (11)

Both the sign from H_normal=-y/tau and the proper derivative scaling
cancel under absolute values. Thus (11) is exactly the norm used in
A.18, not a conformal-coordinate derivative estimate. At the largest
allowed delta, lambda=31/18000000000, the errors are

    (31/28125000000,93/781250000,
     3131/292968750,149048/146484375),

strictly below (1/100,1/2,3,16). Their differences from those caps are
positive exact fractions, replayed independently. Every smaller delta
improves them linearly. The actual metric/state pair therefore belongs
to the anchored radiation slice of the pinned history tube. No claim
that all reference powers p in [1/2,2/3] are realized by this pure
thermal photon state is made.

## 6. Regular open spacetime and the finite mathematical endpoint

As y decreases to zero, Phi(y) tends to positive infinity. Equation
(8) therefore covers x down to negative infinity. At the other end,
y approaches ycrit from below, with finite Phi(ycrit). Its endpoint is

    x_end=integral_2^ycrit (1-2lambda y^2)/
                                  [2y^2(1-lambda y^2)] dy.

The integrand is positive below ycrit and strictly smaller than
1/(2y^2), so 0<x_end<1/4. The named lambda is less than 1/64, so
ycrit>4. For y in [2,4], D>1/2 and 1-lambda y^2<=1, and hence the
integrand is greater than 1/(4y^2). Its integral from 2 to 4 is 1/16.
This proves the useful exact interval

    1/16<x_end<1/4.                                      (12)

For every fixed positive delta, (7) has the finite positive endpoint
value a_end^4=16lambda(1-4lambda). The proper scale factor is smooth
and positive on the open product spacetime x in (-infinity,x_end)
times R^3. Its conformal clock is smooth and strictly monotone there,
so the state construction and the local Hadamard/SEE equations remain
valid at every interior point. The singular endpoint is not part of
the field/metric domain. As delta tends to zero its limiting positive
scale tends to zero; no uniform positive endpoint scale is claimed.

The denominator D tends to zero, while the numerator of F in (6)
has a nonzero positive limit. Thus y' tends to positive infinity and
Hdot_normal tends to negative infinity. Direct curvature calculation gives

    tau^2 R_FK=6(-F+2y^2)
                =-12lambda y^4/(1-2lambda y^2)
                  ->-infinity.                           (13)

There is no smooth or C2 continuation through this regular-coordinate
endpoint with the same metric; the invariant curvature divergence is
not a removable coordinate singularity. The inextendible comoving
normal in the stated low-branch spacetime has finite remaining proper
length tau*x_end. Hence this mathematical branch is timelike
geodesically incomplete in the contracting direction (past-directed
for the expanding interpretation).

This endpoint is not a new QEI-only proof: from the actual SEE,
kappa tau^2 EED=3(F-y^2)>0 on the low branch, so this state already
satisfies timelike convergence. Also H^2 reaches 1/(2b), an
anomaly/gravitational scale; no control of fundamental quantum gravity
or neglected gravitational corrections is asserted. The conditional
future caps in A.18 are not evaluated through this endpoint or claimed
on every shorter incomplete segment. A.18 is an independent
all-Hadamard theorem; the present actual thermal history is an optional
source-compatible realization, not an extra premise of its proof.
