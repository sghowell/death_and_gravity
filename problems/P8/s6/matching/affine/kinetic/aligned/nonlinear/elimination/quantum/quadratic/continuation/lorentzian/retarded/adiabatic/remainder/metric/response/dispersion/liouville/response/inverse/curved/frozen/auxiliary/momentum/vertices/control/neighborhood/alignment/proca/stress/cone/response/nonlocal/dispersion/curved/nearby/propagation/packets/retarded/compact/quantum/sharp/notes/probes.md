# Complete L2 remainder and new rational-width compact probes

Let R(t,s,x)=G-G_clock_front-G_matter_front for ordered endpoints
in the inner interval. The fronts and their positive clock lower
bound are unchanged. Use convolution Fourier normalization
(2*pi)^-3, hence squared L2 norm has Plancherel factor (2*pi)^-3.

Split every radial momentum into exactly three regions.
On 0<=k<=1, abs(Ghat)<4 and the sum of both front multipliers
is at most 2*3*4T=24T, since each integrated radius is <=4T.
Thus abs(Rhat)<5 and this region contributes at most
25/(6*pi^2) to the squared spatial L2 norm.

On 1<=k<=K=10^14, abs(Ghat)<100/k and each front amplitude is
below three, so abs(Rhat)<110/k. Its squared contribution is
at most 110^2*(K-1)/(2*pi^2).
On k>=K, the new absolute error C/k^2, C=10^16, gives
C^2/(2*pi^2*K). No frequency range, low-frequency instability
bound, source factor or phase zero is omitted.

Using pi>3 and enlarging K-1 to K gives the rational upper
25/54+110^2*K/18+C^2/(18K)<10^18.
Therefore the complete remainder has L2 norm below Rbar=10^9,
uniformly in both times. Its continuity and distributional
meaning follow as in S6.91, with these smaller dominating bounds.

Define a NEW pair J_sharp,f_sharp. Retain the S6.92 bump
beta(z)=exp(1-1/(1-z^2)) on |z|<1, zero otherwise, and both
actual integrals B1=int beta and B3=int_R3 beta(|y|).
Keep T=10^-7, s0=-T/4, t0=T/4, h=99/(8*10^16), a=h/10,
Smin=T/16, Smax=3T, Aclock_min=10^-8 and
D=Aclock_min*Smin*a^3/(1024*Smax^2).
Put E0=(D/(32a^2))^2 and the NEW widths
epsilon=E0/10^18, rho=epsilon/20.
These are positive rationals, not the old symbolic-exponential
widths. The code instantiates both literal piecewise bump
formulas at these widths and the fixed source/detector centers.

Explicitly J=rho^-4 beta((s-s0)/rho) beta(|y|/rho)/(B1 B3).
The detector is the product of four bumps with arguments
(t-t0)/a, |x_perp|/a, (x1-S_c(t,s0))/a,
and (|x|-S_c(t,s0))/epsilon, with
S_c(t,s0)=integral_{s0}^t omega_c(u)du.
It is defined by zero outside its interior time window.

Every S6.92 support, cap, plateau and Jacobian inequality is
rechecked with the new epsilon and rho. In particular epsilon<a/100,
all detector points remain within 6a<h of the original detector
center, and the source lies within rho<h of its original center.
Every source/detector support pair remains physically matter-spacelike
and strictly time ordered. The radial detector is away from the
origin. The same all-orders bump proof gives real smooth compact
tests; B1 and B3 are not replaced by one.

The positive clock-front cap still contributes at least D;
the matter shell contributes zero. The norm estimates are
norm(J)_{L1_s L2_y}<=8 rho^-3/2,
norm(f)_{L1_t L2_x}<=8 a^2 sqrt(epsilon).
Because J is nonnegative with spacetime L1 norm one, translating
the spatial L2 remainder and then integrating the source gives
abs(<f,RJ>)<=Rbar*8a^2 sqrt(epsilon)=D/4.
Thus the new full classical pairing is at least 3D/4.

For source-only spatial projection to |k|<=Lambda, use the
all-k>=1 bound abs(Ghat)<=100/k. Plancherel and the full source
and detector norms bound the omitted pairing by
6400*100*a^2/(Lambda*epsilon).
Choose the new exact rational cutoff
Lambda=1+25600*100*a^2/(D*epsilon)<10^94.
The tail is strictly below D/4 and the retained pairing is
at least D/2. The projected source is NOT spatially compact.

The new unprojected smooth compact tests are valid in S6.93's
unchanged CCR algebra. Its source-volume convention therefore gives
[O(f_sharp),O(V J_sharp)]=-i*hbar*C_sharp/kappa times identity,
where C_sharp>=3D/4. This retains the same state-independent
commutator lower, now for this explicitly different pair.
