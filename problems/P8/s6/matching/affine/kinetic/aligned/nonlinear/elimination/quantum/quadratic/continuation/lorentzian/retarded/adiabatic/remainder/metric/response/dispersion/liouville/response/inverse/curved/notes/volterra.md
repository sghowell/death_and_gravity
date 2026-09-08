# Prepared coupled inverse by a weighted Volterra equation

All equations are linearized about the fixed background and
state. This is not a nonlinear existence theorem. Keep the
two-current normal form proved in matching.md, and denote
its full retained Gaussian-plus-fixed-tadpole operator by Q.

## Eliminate the classical matter charge without dividing by lapse

The literal S6.57 density, divided by a^3, is

L=-3v'^2+(J_e+w^2/2-3Theta^2)n^2
  +6Theta*n*v'+s'^2/2+w*n*s'-3ell*v'*s.

The identity ell'=-3H*ell implies that the matter Euler
equation is equivalent to

[a^3*(s'+w*n+3ell*v)]'=0.

Prepared zero data set this conserved quantity to zero.
Thus s'=-3ell*v-w*n. Inserting this into the two metric
Euler equations gives the local system

E_n=2*(J_e-3Theta^2)*n+6Theta*v'-3ell*w*v,
E_v=6v''+18H*v'-6Theta*n'
    -(6Theta'+18H*Theta+3ell*w)*n-9ell^2*v.

Equivalently these are the Euler equations of

L_eff=-3v'^2+6Theta*n*v'+(J_e-3Theta^2)n^2
      -3ell*w*n*v-(9/2)*ell^2*v^2.

The code independently derives the charge, both metric
equations and this effective density from the literal frozen
action. There is no inverse Theta,H,q or J_e in this step.
The matter displacement is recovered by one zero-past
integral afterward.

Pass to physical f=(n,zeta) with v=zeta-delta*n,
delta=1/(2h). If T=[[1,0],[delta,1]], the physical force
operator is T^(-T) E_regular T^-1, with all time-dependent
coefficients kept inside derivatives. This has order at most
two; its highest matrix is
6*[[delta^2,-delta],[-delta,1]]. Both direct physical Euler
variations are checked. No classical principal inverse is
needed.

The classical background equations already vanish. The
vector background gradient is cancelled by the fixed scalar
tadpole. Hence the full physical force is (-delta rho,
3 delta p), with the S6.69 nonlinear-chart contact cancelled
only in the combined response. The quantum operator here
includes that fixed tadpole Hessian, not a reselected profile.

## The equation and the actual preconditioner inverse

In the convention of FORMULATION.md the two metric equations
with a prepared external two-force are

E f + gamma Q f = g, gamma=1/(64*pi^2*L^2)>0.

Four primitives and division by gamma give

(P+V+gamma^-1 I4 E)f=gamma^-1 I4 g.

The local order of E is at most two, so W=V+gamma^-1 I4 E
is a Volterra integral operator with a uniform scalar
majorant v(r)=C*(1+abs(log r)), 0<r<=1, for some finite C.
Large gamma^-1 enlarges C; it does not change integrability.

S6.72 gives the exact inverse Rpre=S^-1 Rop_1 S^-1 of P.
It is an instantaneous smooth matrix A(t) plus a causal
kernel U(t,s) bounded by C_S*||Kreg(t-s)||. The latter
majorant is in L1(0,1). The instantaneous part is nonzero
and is kept. Thus the equation is exactly

f+K f=b,
K=Rpre W, b=gamma^-1 Rpre I4 g.

K is still a Volterra integral operator: its kernel is

A(t)W(t,s)+integral_s^t U(t,r)W(r,s)dr.

It is bounded by an integrable scalar majorant

w=||A||_infinity*v+C_S*(||Kreg||*v)

on (0,1). Young's elementary integral inequality proves its
integrability. This step does not assume the full Rpre norm
vanishes on short intervals; only the genuinely integral W
needs the small-lag property.

## Full fixed interval, not just a formal iteration

Use the norm
||f||_lambda=sup_(t in I) exp[-lambda*(t-u0)] ||f(t)||.
Then

||K||_lambda <= integral_0^1 exp(-lambda*r)*w(r)dr.

Since w is integrable and has no delta term, dominated
convergence makes this bound tend to zero as lambda tends
to infinity. Choose a finite lambda so it is below 1/2.
The series sum_(j>=0)(-K)^j converges in that Banach norm,
giving a unique continuous solution on the entire I.
The weighted norm is equivalent to the ordinary uniform
norm on I; therefore the inverse has a finite C0 bound.

For example a valid symbolic bound is

||f||_C0 <= 2*exp(lambda)*gamma^-1
 ||Rpre||_(C0 to C0) * ||I4 g||_C0.

Also ||I4 g||_C0<=||g||_C0/24 on the unit-length interval.
The proof supplies finite lambda and inverse norm, not
their numerical values. It provides neither a small
one-loop correction nor an optimized continuation window.

Every term is causal. If g vanishes before some preparation
time, uniqueness on that initial subinterval makes f vanish
there too. The Neumann solution is not obtained by dropping
growing pole residues from a formal full propagator.

## Smoothness and the original equations

All fixed diagonal derivatives of W retain an integrable
majorant by matching.md. The same is true for U, since
(partial_t+partial_s) annihilates its lag-only Kreg factor
and only differentiates the smooth multipliers. Under the
causal kernel composition above, translate the integration
variable with both endpoints; the resulting derivatives
have convolution majorants in L1.

For a prepared f, differentiating a Volterra integral can
therefore be written with the derivative on f plus the
diagonal derivative of its kernel. There are no lower
endpoint terms because all initial jets vanish. Iterating
gives

(I+K) f^(n)
 = b^(n)-sum_(j=1)^n binomial(n,j) K_[j] f^(n-j),

where K_[j] has the jth diagonal-derivative kernel. A
difference-quotient argument with the integrable majorants
justifies the first derivative; induction gives every
higher derivative using the same bounded inverse I+K.
For smooth prepared g, b is smooth and prepared because
Rpre's lag convolution can be differentiated onto its
smooth prepared input. Thus f is smooth and prepared.

Finally differentiate the integrated identity four times
in the zero-past distribution algebra. It recovers exactly
E f+gamma Q f=g. For smooth prepared f, S6.68's original
finite mode-response definition is available, so this is
also the original smooth retained response equation, not
only an auxiliary integral equation. Conversely, any smooth
prepared solution of those equations satisfies the
integrated Volterra equation and is therefore the same
solution. The reconstructed s solves the original free
matter equation and its zero-charge condition.

If external forces are represented covariantly, their clock
component is fixed by the linearized Ward identity as in
S6.57; it is not another independently prescribed source.
No claim about an arbitrary unconstrained spacetime tensor
forcing or arbitrary initial-state data is made.
