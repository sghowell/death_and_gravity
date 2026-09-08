# The retarded kernel really is an L1 convolution kernel

This argument supplies the function-space step not obtained
from real-axis matrix bounds alone. All constants below can
depend on fixed positive h and m. Compact positive parameter
sets are uniform by the explicit formulas and gaps. No
numerical bound on the full L1 norm is asserted.

## Differentiable endpoint estimates

Put a(Omega)=rho(Omega^2), Omega>=Omega0=2m. At threshold the
closed cut formula has

a(Omega)=Cth*sqrt(z)+O(z^(3/2)), z=1-Omega0^2/Omega^2,
Cth=T^-1 M(0) T^-1/8.

More precisely rho(z)=sqrt(z)*A(z), where A is real analytic
near zero: sqrt(z)*atanh(sqrt(z)) is analytic in z, the density
numerator is sqrt(z) times an analytic matrix, and its
denominator is nonzero at zero. Thus, with x=Omega-Omega0,

a(Omega0+x)=A0*sqrt(x)+O(x^(3/2)),

with the same expansion differentiable termwise. A0=Cth/sqrt(m).
The certificate records Cth exactly; it is nonzero rank one.

At infinity let l=log(Omega^2/m^2). Schur inversion on the cut
gives

a(Omega)=Cinf/l^2+O(l^-3)+O(Omega^-2*polylog(Omega)),
Cinf=g*g^T/2,
g=Q*(-b*/a*,1)^T=(-12231h/3904,395/244)^T.

The exact expressions are rational in sqrt(1-e), e and
atanh(sqrt(1-e)), e=4m^2/Omega^2. The last function is
log(Omega/m) plus a convergent analytic power series in e.
The nonzero a* and logarithmically growing Schur denominator
justify differentiating the inverse expansion. In particular,

||a(Omega)||=O(log(Omega)^-2),
||a'(Omega)||=O(1/[Omega*log(Omega)^3]),
||a''(Omega)||=O(1/[Omega^2*log(Omega)^3]).

These derivative estimates are consequences of the explicit
closed expressions, not differentiation of an unspecified
big-O remainder. The prime here differentiates Omega.
In particular a is bounded, a' is absolutely integrable on
[Omega0,infinity), a(Omega0)=a(infinity)=0.

## Oscillatory integral and small time

Define, for t>0,

Kreg(t)=-2 integral_(Omega0)^infinity a(Omega)*sin(Omega*t) dOmega.

This improper integral exists by one integration by parts
using a' in L1 and a(infinity)=0. It is continuous for t>0,
and equals its Abel limit with e^(-epsilon*Omega) inserted.
The integral without the sine is divergent at infinity, so
absolute frequency convergence is expressly not used.

For small t take R=1/t sufficiently large. The tail is bounded
by

(1/t)[||a(R)||+integral_R^infinity ||a'(Omega)||dOmega]
 =O(1/[t*log(1/t)^2]).

For the initial portion use abs(sin(Omega*t))<=Omega*t.
Below sqrt(R), boundedness of a gives O(1); above sqrt(R),
||a||<=C/log(R)^2 gives O(t*R^2/log(R)^2).
Hence the same integrable bound holds for the whole kernel:

||Kreg(t)||=O(1/[t*log(1/t)^2]).

The integral of this bound at zero is finite, as seen by
v=log(1/t). An O(1/[t*log(1/t)]) bound would not suffice.

## Large time and the threshold tail

Choose a smooth compactly supported cutoff chi(x) equal to
one near x=0. Subtract A0*sqrt(x)*chi(x) from a(Omega0+x).
The remainder b is C1 at zero, b(0)=b'(0)=0 and b'' is in L1:
near zero its second derivative is O(x^-1/2), on compact
positive x it is smooth, and at infinity use the above
a'' estimate. Also b and b' tend to zero at infinity.
Two integrations by parts therefore bound its Fourier
integral by ||b''||_L1/t^2.

For sqrt(x)*chi(x), split at x=1/t. Its initial integral is
O(t^-3/2); twice integrating the tail by parts gives endpoint
terms O(t^-3/2) and a second-derivative integral
O(t^-2*integral_(1/t)^1 x^-3/2 dx)=O(t^-3/2).
The extra phase exp(i*Omega0*t) has unit modulus.
Consequently Kreg(t)=O(t^-3/2) as t tends to infinity.

Together, the two endpoint estimates prove
Kreg in L1((0,infinity), real symmetric 2x2 matrices).

## Laplace transform without an invalid Fubini step

For epsilon>0 define a_epsilon=e^(-epsilon*Omega)*a.
Its frequency integral is absolutely convergent. Fubini
with the exponentially damped time integral gives

Laplace[Kreg_epsilon](s)
 =-integral_(Omega0)^infinity
  2*Omega*e^(-epsilon*Omega)*a(Omega)/(s^2+Omega^2) dOmega,

initially for real s>0. This spectral integral has an
epsilon-independent integrable majorant at infinity
O(1/[Omega*log(Omega)^2]).

The preceding time-domain bounds can also be chosen uniformly
for 0<epsilon<=1. For the small-time tail, the added derivative
term has integral at most sup_(Omega>=R)||a(Omega)||, since
integral_R^infinity epsilon*exp(-epsilon*Omega)dOmega<=1.
For large time, subtract the damped threshold leading term.
The remainder's second derivative is uniformly integrable:
besides the existing second derivative, its extra terms
are bounded by epsilon*||a'||_L1 and
epsilon^2*integral exp(-epsilon*Omega)||a||dOmega,
uniformly finite. On the fixed threshold neighborhood the
analytic remainder bounds are uniform in epsilon. Thus
dominated convergence is legitimate in time as well.

Taking epsilon down to zero and changing tau=Omega^2 yields

Laplace[Kreg](s)=-integral_(4m^2)^infinity rho(tau)/(s^2+tau) dtau.

The formula extends to Re(s)>0 by holomorphy. Together with
the previous note this proves

Laplace[Rinf*delta+Kreg](s)=B(s^2)^-1.

B(s^2) has logarithmic growth and is itself the Laplace
transform of a zero-past distribution; equivalently it is
defined by the exact subtracted retarded bubble plus F.
Multiplication of these Laplace transforms proves both
convolution inverse identities in the zero-past distribution
algebra. No extra initial homogeneous solution or pole
residue has been removed. On continuous inputs, convolution
with Kreg is an ordinary L1 convolution.

## Domains, moment bound and parameter scaling

For bounded continuous f on a half-line, or C0([0,T]),

Rop f(t)=Rinf*f(t)+integral_0^t Kreg(t-v)f(v)dv

is continuous and bounded in the same norm by
||Rinf||+||Kreg||_L1. Here Rinf*f(t) means matrix multiplication,
not another time integral. Smooth inputs zero on an initial
neighborhood stay so, since derivatives can be put onto f.
A generic continuous input is not claimed to have a C1 output.

The regular-kernel norm on (0,T) tends to zero with T, but
the full operator norm does not: Rinf is nonzero. It is not
valid to apply a scalar small-inverse contraction to both
channels without separating this instantaneous complement.

The absolutely convergent primitive is

J(t)=integral_0^t Kreg(v)dv
 =-integral_(4m^2)^infinity
   rho(tau)/tau*[1-cos(sqrt(tau)*t)]dtau.

The static sum rule gives -2*Mmoment<=J(t)<=0 in Loewner
order, hence ||J(t)||_2<=2*trace(Mmoment).
For f in C1 with f(0)=0, integration by parts gives

Rop f=Rinf*f+J*f',
||Rop f||_C0<=trace(Rinf)||f||_C0
               +2*trace(Mmoment)*T*||f'||_C0.

The certificate provides these exact trace bounds for h=1
and h=125/64. This C1 estimate is additional information,
not the proof of the C0 endomorphism.

At fixed h, B_m(s^2)=B_1((s/m)^2), so
Kreg_m(t)=m*Kreg_1(m*t). Its full half-line L1 norm is therefore
independent of m>0. No uniform massless limit or m=0 inverse
is inferred.
