# Original causal inverse on a specified all-momentum flat TT graph

All derivatives below are causal distributional derivatives. Initial atoms are retained. The force convention is O_phys h=f, not the negative unsourced Euler residual.

For q=P^2>=0 put L_q=D^2+q and s_q(t)=sin(sqrt(q)t)/sqrt(q), with s0(t)=t. The original Bromwich inverse of the complete reciprocal is

G_phys(t,q)=C0 theta(t){s_q(t)/C
 +2 Re[R sinh(sqrt(z-q)t)/sqrt(z-q)]
 +integral rho_E(tau) s_(q+tau)(t) dtau}.                   (1)

The square root in the displayed term has positive real part. Both conjugate first-sheet poles are retained. Since Re sqrt(z-q)<=sqrt(|z|)<sqrt(kappa) for every q>=0, a Bromwich line Re lambda>sqrt(kappa) lies to their right uniformly. No frequency-sign-dependent sheet selection or pole-excluding Fourier prescription is made.

## Uniform bounds, without frequency division

The real-q sine kernel and its derivative have bounds t and1. Duhamel about L_q for the complex kernel gives the absolutely convergent Volterra series with bounds

|s_complex(t)|<=sum_(n>=0) |z|^n t^(2n+1)/(2n+1)!
 <=t exp(sqrt(|z|)t),
|s_complex'(t)|<=sum_(n>=0) |z|^n t^(2n)/(2n)!
 <=exp(sqrt(|z|)t).

This avoids any division by a difference of nearly equal frequencies. The exact moment/residue estimate in notes/dispersion.md and C0/C<1 give

|G_phys(t,q)|<8t exp(sqrt(kappa)t),
|D G_phys(t,q)|<8 exp(sqrt(kappa)t)                          (2)

for t>0, uniformly in ALL spatial q. The kernel vanishes at t0; its first derivative has no extra delta. On any finite window of length T, Minkowski and Young inequalities give both continuous-time and L2-time bounds

||G_phys||<=4T^2 exp(sqrt(kappa)T),
||D G_phys||<=8T exp(sqrt(kappa)T)                           (3)

on every H^r_TT, without spatial derivative loss. These huge finite bounds are not stability.

## The actual forward operator and its weak target

Use the original paired representation, without splitting divergent local/cut pieces:

A2_q=delta/30+L_q H_q,
H_q(t)=theta(t) integral_0^1 [W2(y)/(1-y^2)]
                       s_(q+4m^2/(1-y^2))(t) dy.

The integral is absolutely convergent and

|H_q|<= [1/(2m)] integral W2(y)/sqrt(1-y^2) dy
      =9pi/(128m)                                          (4)

uniformly in q and time. The beta moment is independently replayed. Therefore the full original force operator is

O_phys=[C L_q+L_q^2/30+L_q^3 H_q]/C0.                       (5)

Equation (5) defines a continuous forward map C_t H^r_TT to D'_t H^(r-6)_TT on finite slabs: use time-test derivatives through order6, the uniform convolution bound (4), and spatial powers q through3. The original causal atoms, including the source germ boundary, are included in this distributional definition.

Fix the original nonempty zero initial neighborhood, not a new reset. Let Y_r be C(I;H^r_TT) with that fixed zero germ and define

Dom O_phys={h in Y_r: the full causal distribution O_phys h is in Y_r}.

The graph norm is ||h||_Y+||O_phys h||_Y. This domain is explicit; no larger maximal distribution domain is asserted.

At bounded Fourier support, the original analytic product O_phys(p) C0 E_C(p)=1 gives both convolution identities in causal distributions, by the Laplace transform and associativity. For a finite interval, causal extensions beyond its right endpoint do not change these identities within the interval.

Now let spatial cutoffs tend to infinity. They commute with both flat operators and converge strongly on Y_r (uniformly in time because a continuous path has compact range). Equations (2)-(3) give uniform convergence of inverse images. The weak forward continuity of (5) then gives O_phys G_phys f=f for every f in Y_r. Conversely, for h in the stated graph, apply the bounded-support identity to each cutoff of h and its full forward value, then pass to the same limit to obtain G_phys O_phys h=h.

Thus the inverse is genuinely two-sided on THIS all-momentum flat TT causal graph, with inverse-to-graph norm at most1+4T^2 exp(sqrt(kappa)T). The same construction supplies the analogous L2-time realization if its domain is defined in that space; it is not a claim about arbitrary distributions.

For smooth time-prepared input, time derivatives commute by transferring them to that input with zero original boundary terms. The solution is smooth in time in the same spatial H^r. No spatial momentum differentiation is needed or claimed.

This is not the curved scalar/clock/matter operator, not S222's unrestricted completed graph, not small physical feedback, and not a stable or nonlinear UV-parent solution.
