# S6.87: actual ordinary-Proca prepared homogeneous coupled inverse

Keep the source-free constant-physical-mass Proca candidate S6.81, the
already-fixed S6.82 profile and the original S6.55 state. The clock,
strict classical margin, fixed m=1000 and L=10^400 are unchanged.
The prepared source interval is I=[-1/2,1/2]. Independent initial-state
variations and arbitrary higher-derivative initial data are excluded.

The S6.86 scalar inverse is not an inverse of the two-source loop block.
This checkpoint instead uses the actual classical constraint direction.
For physical lapse/log-scale sources (n,v), define the prepared bijection

    n=eta', v=w+H eta, eta(initial)=0.

The regular clock scale is vhat=w+H eta-delta eta', delta=1/(2h).
Both changes are regular at H=Theta=0. The free matter charge is set to
its unchanged prepared value, not divided by H, Theta or spatial momentum.

Let T be the actual matter-reduced, normalized tree operator in (eta,w).
Let Q be 64*pi^2*L^2 times the full retained vector-plus-fixed-profile
force response, transformed with E_eta=-(partial_u+3H)E_N+H E_Z,
E_w=E_Z. Put gamma=1/(64*pi^2*L^2). All coefficient multipliers remain
inside their time derivatives. For original prepared force g=(g_N,g_Z),
set g_eta=-(partial_u+3H)g_N+H g_Z and g_w=g_Z.

## Normal form and inverse

Exact mode time covariance, all adiabatic orders and the finite local
stress show that the time channel of Q is local of order at most two.
Its only nonlocal channel is w. The actual tree has derivative-order
matrix [[4,3],[3,2]], with fourth time coefficient A(u)=-6 delta(u)^2.
On the entire interval, |A|>=6144/15625>0.

For four zero-past primitives I4 and the exact new scalar block F_m,

    I4(T+gamma Q)=diag(A(u),gamma F_m(partial_u^2))+V.

V is an actual causal integral operator with a weak-log integrable
kernel; all fixed diagonal derivatives obey the same class. Its first
row is stronger: a bounded local-primitive kernel with bounded first
output derivative. The new leading dimensional current, finite local
coefficient and all-momentum subtraction are matched explicitly.
No old full-rank mass-profile inverse or arbitrary local fourth contact
is transferred.

The diagonal inverse is multiplication by A^-1 and convolution by
gamma^-1 times the new L1 scalar kernel. Composition with V gives an
integrable Volterra majorant. A finite exponentially weighted norm
therefore makes the genuine integral remainder contractive on all I.

Every smooth prepared original two-force has a unique smooth prepared
solution of the retained homogeneous tree-plus-Gaussian-plus-fixed-profile
equations. It is causal and obeys a finite C0 bound in the original
physical sources; the first-row local structure controls n=eta' without
requiring a derivative of arbitrary C0 w. No numerical inverse norm,
small correction, stability or quantum-cone estimate is asserted.
The matter displacement is recovered with the original negative-sign
zero-charge reconstruction. The force transformation adds no independent
prepared solution because its remaining homogeneous kernel is a
zero-initial first-order equation.

## Boundary

This is retained linear response, not a nonlinear quantum bounce, an
arbitrary spatial/initial-state theorem, a physical all-frequency EFT
propagator, Wilson matching, UV completion or original P8 closure.
A causal finite-interval inverse can have a very large norm and does
not exclude growing modes or establish healthy quantum propagation.

Evidence is exact algebra, unchanged-state high-frequency estimates,
common dimensional matching and written Volterra proofs. The analytical
existence arguments are not proof-assistant formalization.
