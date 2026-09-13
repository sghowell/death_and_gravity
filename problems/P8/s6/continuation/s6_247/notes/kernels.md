# Complete two-threshold causal kernels

Let a_i(Omega)=rho_i(Omega²), Omega0=2mu, OmegaH=2sqrt(n). The full density is positive between and above these thresholds, vanishes at Omega0 and at infinity, and is continuous at OmegaH. Its interior value at OmegaH is not zero.

The closed formulas and the lowest-threshold gap give a_i(Omega0+x)=A_i sqrt(x)+O(x^(3/2)) with a convergent local expansion. Near OmegaH the trace density has BOTH left and right square-root terms derived in notes/measure.md, after the analytic change from log spectral variable to frequency. The shear's heavy-threshold nonanalytic term begins at order|Omega-OmegaH|^(5/2). On compact intervals away from these two points the complete density is smooth.

At high frequency use u=log(s/(4m²)) separately for each mass and the stable exact identity
atanh(sqrt(1-exp(-u)))=u/2+log(1+sqrt(1-exp(-u))).
The full polynomials and their derivatives give
a_i=O((log Omega)^-2),
a_i'=O(1/[Omega(log Omega)^3]),
a_i''=O(1/[Omega²(log Omega)^3]).
All constants may depend on the actual fixed masses. This is obtained by differentiating explicit smooth formulas, not an unproved remainder symbol.

It follows that a_i' is integrable on[Omega0,infinity), including both sides of the interior threshold. The oscillatory expression

K_i(t)=-2 theta(t) int_(Omega0)^infinity a_i(Omega)sin(Omega t)dOmega

is therefore defined by one integration by parts for every t>0, uniformly on each compact positive-time interval. The interior integration boundaries cancel because a_i itself is continuous. No unweighted absolutely convergent frequency integral is asserted.

## Integrability near zero and infinity in time

For sufficiently small t set R=1/t so sqrt(R) exceeds both physical thresholds and the fixed start of the logarithmic tail. Split the frequency integral at sqrt(R) and R.

In the low region, boundedness of a and |sin(Omega t)|<=Omega t gives O(1). In the middle, a<=C/log²R gives O(1/[t log²R]). A single tail integration by parts gives
[|a(R)|+int_R^infinity |a'|]/t
with the same bound. Fixed earlier intervals only change the constant. Thus K_i is integrable near t0.

For large t, subtract localized one-sided models for the first square-root threshold and BOTH interior trace square roots. At the interior threshold use a convergent one-bank Puiseux expansion of the reciprocal, with the nonzero open-Proca denominator. Integer-power terms have the same smooth continuation on both sides; after subtracting the leading half powers, the residual has continuous first derivative and an integrable second derivative. The shear interior term already has integrable second derivative.

A smooth cutoff equal to one near each threshold localizes its model. Compare each one-sided sqrt(x) model with sqrt(x)exp(-x); the difference is smooth away from x0 and has integrable second derivative at x0. The latter model has exact Fourier transform Gamma(3/2)(1-/+it)^(-3/2), with the appropriate oscillatory threshold phase. Each threshold contributes O(t^-3/2). The remaining full density has two integrable derivatives and contributes O(t^-2).

Hence K_i belongs to L1(0,infinity), for BOTH total channels. No numerical value of these L1 norms is claimed. This differs from the old Proca-only shear reciprocal, whose isolated undamped pole prevented half-line L1. The changed result follows from the verified zero set of the specified sum, not a deletion of that old pole.

## An absolutely convergent shifted primitive

For every external q>=0, define

J_i,q(t)=-int_(4mu²)^infinity rho_i(s)/(q+s)
                    [1-cos(sqrt(q+s)t)]ds.

The positive shifted static measure has mass1/A_i(q)<=1/c_i. This gives continuous J_i,q, J_i,q(0)=0, and
-2/c_i<=J_i,q<=0.
In particular
|J_trace,q|<=1/(ell+2), |J_2,q|<=120/(ell+2).

The primitive integral is absolutely convergent, including the full high-frequency tail. Its Laplace transform follows by Fubini against this finite static measure:

L[J_i,q](lambda)=1/[lambda F_i(lambda²+q)], Re lambda>0.

For q0, its derivative is the above ordinary K_i, also as a causal distribution at the initial boundary. The all-q matrix construction uses J directly and needs no unproved uniform scalar-kernel L1 norm.

## Paired full forward distributions

For each species of mass m, let c=1-y² and Omega_y,q=sqrt(q+4m²/c). The full forward kernel is

G_species,i,q(t)=theta(t) int_0^1
 W_species,i(y) sin(Omega_y,q t)/(c Omega_y,q)dy.

The integrable pointwise majorant W/(2m sqrt(c)) gives bounded continuous G with G(0)=0. Exact beta moments give the heavy bounds41pi/(64sqrt(n)) for trace and pi/(384sqrt(n)) for shear; the old Proca bounds remain27pi/(64mu) and9pi/(128mu).

Define
F_i,q=-c_i delta-(D²+q)(G_Proca,i,q+G_heavy,i,q).
Its Laplace transform is exactly F_i(lambda²+q), with both full mass integrals. This paired expression does not separate two individually divergent local and spectral pieces. All initial distributional terms are retained.

At q0, transform uniqueness yields F_i*K_i=K_i*F_i=delta. Thus each ordinary factor inverse is valid on its full zero-past causal graph, not merely as an equation at t>0. The ordinary kernel and all fixed-window graph identities use only the specified flat reference; they are not yet an actual curved constrained inverse.
