# The scalar inverse has an ordinary L1 causal kernel

This applies the generic oscillatory-kernel argument proved in the [S6.72 kernel note](../../../../../../../../../../../../../../../notes/causal-kernel.md), not that checkpoint's different matrix inverse. All needed endpoint hypotheses are re-established here from the new scalar density. Constants may depend on fixed m>0; no numerical full L1 norm is asserted.

Put a(Omega)=rho(Omega²), Omega>=2m. The explicit cut formula gives a(2m+x)=A0 sqrt(x)+O(x^(3/2)), A0=675/(512 sqrt(m)), with a differentiable analytic remainder. At infinity its rational/logarithmic form gives

    a=O(log(Omega)^-2),
    a'=O(1/[Omega log(Omega)^3]),
    a''=O(1/[Omega² log(Omega)^3]).

These follow by differentiating the explicit formula, not an unqualified big-O term. In particular a is bounded, a' is integrable and both endpoint values of a are zero.

For t>0 define K(t)=-2 integral_(2m)^infinity a(Omega) sin(Omega t) dOmega. One integration by parts using a' in L1 gives an improper, continuous oscillatory integral. The frequency integral without its oscillation is not absolutely convergent.

For small t take R=1/t. The high-frequency tail is bounded by a constant times [|a(R)|+integral_R^infinity |a'|]/t, hence O(1/[t log(1/t)^2]). Below sqrt(R), use boundedness and |sin(Omega t)|<=Omega t to obtain O(1). Between sqrt(R) and R, use the logarithmic bound to obtain the same O(1/[t log(1/t)^2]) estimate. This time singularity is integrable at zero.

For large t subtract A0 sqrt(x) times a smooth compact cutoff equal to one near x=0. The remainder and its first derivative vanish at zero, tend to zero at infinity, and its second derivative is integrable. Two integrations by parts give an O(t^-2) remainder. Splitting the cutoff sqrt(x) term at x=1/t and integrating its tail twice gives O(t^-3/2). Thus K=O(t^-3/2) at infinity, and K belongs to L1(0,infinity).

To justify its Laplace transform, first insert exp(-epsilon Omega) in the frequency integral. Fubini is then valid. The spectral transform is dominated as epsilon decreases to zero by the integrable 1/[Omega log(Omega)^2] tail. The preceding time bounds can be chosen uniformly: the added first-derivative term has tail at most sup_(Omega>=R)|a(Omega)|; the damped threshold subtraction has a uniformly integrable second-derivative remainder. Dominated convergence in time and frequency gives

    Laplace[K](s)=-integral_(4m²)^infinity rho(tau)/(s²+tau) dtau=1/F(s²),

first for real s>0, then throughout Re(s)>0 by holomorphy. The exact subtracted scalar block is itself a zero-past distribution. Multiplication of Laplace transforms proves both scalar convolution inverse identities. No pole residue or independently chosen homogeneous solution is discarded.

For bounded continuous f, convolution K*f is a C0 endomorphism with finite norm bounded by ||K||_L1. Its norm on [0,T] tends to zero as T decreases to zero because there is no instantaneous delta term. A generic C0 input is not claimed to acquire a full derivative; smooth prepared inputs remain smooth and prepared.

The primitive is absolutely represented by

    J(t)=-integral rho(tau)/tau [1-cos(sqrt(tau)t)] dtau.

Its exact static moment gives -1/2<=J(t)<=0. For f(0)=0, integration by parts yields Rop f=J*f' and ||Rop f||_C0<=T ||f'||_C0/2. This additional C1 estimate is not substituted for the L1 proof of the C0 bound.

Finally F_m(p)=F_1(p/m²), hence K_m(t)=m K_1(mt). Its full half-line L1 norm is independent of positive m. No m=0 inverse or massless pointwise limit is inferred. These results concern the scalar normalized range factor; the full loop term still has its fourth-derivative prefactor, and the inactive two-source direction remains outside this inverse.

