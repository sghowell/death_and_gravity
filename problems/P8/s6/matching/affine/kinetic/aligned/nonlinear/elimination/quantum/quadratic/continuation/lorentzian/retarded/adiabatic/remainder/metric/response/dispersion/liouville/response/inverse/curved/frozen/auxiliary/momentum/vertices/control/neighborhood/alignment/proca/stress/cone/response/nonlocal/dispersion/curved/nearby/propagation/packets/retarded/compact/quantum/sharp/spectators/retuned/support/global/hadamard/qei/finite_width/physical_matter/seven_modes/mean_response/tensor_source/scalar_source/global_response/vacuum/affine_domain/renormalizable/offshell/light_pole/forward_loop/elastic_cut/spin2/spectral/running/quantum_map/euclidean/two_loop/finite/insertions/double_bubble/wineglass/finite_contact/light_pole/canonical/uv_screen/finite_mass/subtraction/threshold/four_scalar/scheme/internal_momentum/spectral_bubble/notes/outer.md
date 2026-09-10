# The finite outer integral and its complex-domain bound

Let u denote spectral mass squared, z the outer-channel invariant,
A=x(1-x), and Delta(z)=x u+1-x-A z. Define the finite bubble difference

    B(z;u)-B(2;u) = -integral_0^1 log[Delta(z)/Delta(2)] dx.

Start in real Euclidean kinematics with the common translation-invariant
regulator and the already paired inner on-shell subgraph. The usual
two-propagator parameterization involves only a real momentum translation
there. Form the displayed finite difference before removing the regulator.
Its parameter formula analytically continues to the domain below; no
unjustified complex translation of a divergent momentum integral is used.

For Re(z)<=4 and u>=T>=16,

    Re Delta >= x(u-4)+(1-x)>0,
    Re Delta-x(u-4)-(1-x)=A(4-Re z)+4x^2.

The logarithmic difference is analytic, its derivative is the integral of
A/Delta, and its second derivative is the integral of A^2/Delta^2. Hence

    |B'(z;u)| <= 1/[2(u-4)],
    |B''(z;u)| <= 1/[3(u-4)^2].

The parameter endpoints follow by continuity; no division by x is needed
at x=0. This domain contains |z-2|<=2, including the required forward
unit disc and the zero-transfer channel.

After summing the two possible inner light-line positions the channel is

    A_channel(z)-A_channel(2)
       = (L^2/Q) integral_T^infinity w(u)[B(z;u)-B(2;u)] du.

Only the constant is subtracted. For u>=16, w(u)<=4C/u and u-4>=u/2.
The derivative bound implies an absolutely integrable majorant 4C/u^2
times L^2 |z-2|/Q. Therefore the finite channel difference obeys

    |A_channel(z)-A_channel(2)| <= 4 C L^2 |z-2|/(Q T).

The twice-differentiated majorant is 16C/(3u^3) times L^2/Q.
These bounds justify spectral integration, differentiation and analytic
continuation of the finite parameter representation. The unsubtracted
outer spectral integral is not assigned a finite value.
