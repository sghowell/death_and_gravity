# Ordinary causal kernel and both finite-window inverse identities

All time variables in the estimates can first be scaled by the fixed positive mass. The isolated inverse is

K_m(t)=-theta(t)[R sin(omega_* t)/omega_*+2 integral_(2m)^infinity a(Omega)sin(Omega t)dOmega],
omega_*=m sqrt(r), a(Omega)=rho(Omega^2).

The continuum frequency integral is ordinary oscillatory, not absolutely convergent. Because a vanishes at both endpoints and a' is integrable, one integration by parts defines it for everyt>0 and uniformly on compact positive-time intervals.

For small dimensionless t set Rcut=1/t. Split the integral at sqrt(Rcut) and Rcut. In the low region use bounded a and |sin(Omega t)|<=Omega t, giving O(1). In the middle region a<=C/log^2(Rcut), giving O(1/[t log^2(1/t)]). In the tail integrate once: [|a(Rcut)|+integral_(Rcut)^infinity |a'|]/t has the same bound. A fixed compact initial frequency interval only changes the constant. Thus the continuum kernel is integrable near time zero.

For large t, write x=Omega-2m and subtract the explicit model

A sqrt(x)exp(-x), A=50625/(59168 sqrt(m)).

The remainder b(x) is O(x^(3/2)) at0, with b0=b'0=0 and b'' locally integrable there. The cut derivative estimates give b'' inL1 at infinity as well; b and b' vanish there. Two integrations by parts therefore bound the remainder transform by O(t^-2). The model has exact complex transform A Gamma(3/2)(1-it)^(-3/2) times exp(2imt), giving O(t^-3/2). This proves that the continuum piece belongs toL1(0,infinity). The change in dimensionful exponential scale only affects the harmless constants.

The pole part is bounded on finite windows and has nonzero absolute integral4R/omega_*^2 per period. Hence the full K_m isL1(0,T) for all finiteT, but is NOTL1(0,infinity): the integral of its absolute value is bounded below by that of the pole minus the finite continuumL1 norm. This is why the old trace half-line theorem does not transfer.

## An absolutely convergent primitive

The total static moment gives an everywhere absolutely convergent primitive

J_m(t)=-R/omega_*^2(1-cos(omega_*t))
       -integral rho(tau)/tau (1-cos(sqrt(tau)t))dtau.

It is continuous, J_m(0)=0 and-60<=J_m<=0, since1-cos lies in[0,2] and the total positive measure has mass30. Truncated spectral primitives converge uniformly; their derivatives converge uniformly on each compact intervalt>0 by the a' bound. Thus J_m'=K_m there. The proven localL1 bound and J_m(0)=0 give the identity also as causal distributions.

For Re(lambda)>0, Fubini is valid for the primitive using the finite static measure. Its Laplace transform is

L[J_m](lambda)=-(1/lambda)[R/(lambda^2+omega_*^2)+integral rho(tau)/(lambda^2+tau)dtau]
             =1/[lambda F2(lambda^2)].

Therefore L[K_m]=1/F2(lambda^2). This primitive argument avoids any unproved exchange of the undamped oscillatory frequency integral with a time integral.

## Construct the forward causal distribution without divergent terms

Let c=1-y^2, Omega_y=2m/sqrt(c), and define

G_m(t)=theta(t) integral_0^1 W2(y) sin(Omega_y t)/(c Omega_y)dy.

The integrable pointwise majorant W2/(2m sqrt(c)) has integral9pi/(128m), as follows from the beta moments pi/4,3pi/16,5pi/32. Hence G_m is bounded and continuous, G_m(0)=0, and its Laplace transform for Re(lambda)>0 is

L[G_m](lambda)=integral_0^1 W2(y)/[4m^2+c lambda^2]dy.

For each fixed right-half-plane lambda the y integral is absolutely convergent; the bounded time majorant justifies Fubini. Define the causal distribution

F=-delta/30-partial_t^2 G_m.

It has Laplace transform F2(lambda^2), including all endpoint distribution terms automatically. This definition does not separate the radial integrand into two individually divergent local/cut expressions. On smooth prepared f with f(0)=f'(0)=0 it agrees with-f/30-G_m*f'', and otherwise the initial distribution terms are retained by convolution.

Causal distributions supported in[0,infinity) have a well-defined associative convolution. The specific distributions here have Laplace transforms on Re(lambda)>0. Their transform product is1: lambda^2 cannot hit the negative cut or its isolated real pole while Re(lambda)>0. Uniqueness of the Laplace transform gives

F*K_m=delta, K_m*F=delta.

## Finite-window graph domain and uniqueness

FixT finite and regard sources f inC[0,T] as ordinary causal distributions, including their right-hand value at0. A continuous zero-past response u has u(0)=0. Its graph domain consists precisely of such u for which F*u agrees with an ordinary source theta f as a distribution near the complete initial boundary and throughout[0,T). No initial delta or derivative of a delta is ignored. Local causal convolution makes the condition independent of any extension beyondT.

For each f, u=K_m*f is continuous, u0=0 and satisfies F*u=f. Conversely K_m*(F*u)=u proves uniqueness on that graph domain. Thus the isolated factor has a bounded inverse fromC[0,T] into continuous responses, with operator bound ||K_m||_L1(0,T). No density statement about this graph domain is needed or asserted. The bound tends to zero on shrinking windows by absolute continuity of the localL1 integral; no numerical globalL1 constant is invented.

If f isC1 and f(0)=0, integration by parts gives K_m*f=J_m*f' and

||K_m*f||_C0<=60T||f'||_C0.

The preparation condition is explicit. For f0 nonzero there is an additional J_m(t)f(0) term. Mass scaling is rho_m(tau)=rho_1(tau/m^2), R_m=m^2R_1 and K_m(t)=mK_1(mt). The finite-windowL1 norm is consequently the mass-one norm on[0,mT], not a mass-independent half-line norm of the full kernel.

This is a scalar reference-factor inverse. It does not establish the mapping, normal form, graph invariance or invertibility of the full nonzero-transfer curved quantum-force system.
