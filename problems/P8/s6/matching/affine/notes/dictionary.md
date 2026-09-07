# Exact coefficient matching and solution lifting

## Signature, the target, and the inverse construction

Use the complete frozen CD witness with X_repo=-x and g_source=-g_repo.
The determinant of the metric is unchanged in four dimensions. Constant
overall sign reversal leaves the Levi-Civita connection and covariant
Ricci tensor unchanged; its scalar trace changes sign. Counting inverse
metrics gives signs (+,+,-,-,+) on the five DHOST quadratic operators.
Thus f=-F2_repo, alpha1,2=A1,2, alpha3,4=-A3,4 and alpha5=A5.
The original free matter term becomes -g_source^{ab} chi_a chi_b/2.
Neither the physical causal cone nor which metric couples to matter is
changed. No affine covariant derivative occurs in that first-derivative
matter action; its connection source is identically zero for arbitrary chi.

For reference, the general principal inverse with A1=0 follows by solving
the source equations for f and alpha1. With p=f1, c=f2 and quartic J4,

    c=2(p-f)/x, J4=(p^2-f^2)/(f*x^2), Delta=2p^2/f.

Substitution in alpha3 gives the first-order identity

    p_x/p-f_x/f=x*alpha3/(4f).

These formulas require x,f,p nonzero. They do not by themselves prove
the full connection Hessian has no additional kernel; the separate
literal-action proof does that for the specified CD candidate.

The actual target has h=(1+u^2)^3, w=h-1-x, f=w/(2h), alpha3=1/(hx).
The positive solution p=sqrt(w/h)/2, normalized p(u,-1)=1/2, satisfies
that equation and p^2=f/2. It gives Delta=1 and J4=(1+x)/(2h*x^2).
Substitution yields all five target operators, including alpha4 and
alpha5; the independent audit derives their repository Ia completion
without importing the primary dictionary. More importantly, the full
connection action independently reconstructs these same five operators
from all ten scalar-Hessian components. Published formulas are a
comparison, not the input to that literal Schur elimination.

The closed original tube lies inside the open set
u real, -6/5<x<-4/5. There w/h is strictly between 4/5 and 6/5,
p and f are positive, x is nonzero, and 8p^2-1>3/5. All selected square
roots are the positive real branches. On the smaller closed tube,

    9/20<=f<=11/20, 9/40<=p^2<=11/40, 8p^2-1>=4/5.

The full rank calculation is therefore valid on an open neighborhood,
not only at the rolling background or in a finite time sample.

## The lower-order source formula is derived, not assumed

Literal distortion elimination, with p_u,p_x,c_u,c_x and J3 independent
when extracting the coefficients, gives

    P=J2+3x(p_u-J3*x)^2/Delta,
    Q1=-2f_u+4p(p_u-J3*x)/Delta,
    Q2=2f_u/x-4(p-3x*p_x)(p_u-J3*x)/(x*Delta).

The last bracket differs from the source's printed p-3p_x. All 64
connection Euler equations and the complete reduced density are
checked before comparison with that source. As a separate control,
pure p(phi,x)R_Gamma has reduced correction 3(grad p)^2/(2p), modulo
-3 Box(p). Its cross coefficient is 6p_u*p_x/p. The printed formula
instead adds a factor 1/x to that cross coefficient. Unit positive x
would hide this test; x=-2 with nonzero derivatives does not.

With p^2=f/2 and Delta=1, write Q1=q. Since p_u=f_u/(4p), the required
cubic is J3=-(q+f_u)/(4*p*x). The remaining condition Q2=2q_x becomes

    q_x+Aq=b, A=1/(2x)+3/(4w), b=-3f_u/(2w).

This equation only determines a local coefficient function of the
scalar and its kinetic invariant. It is not an integral over physical
time, a retarded prescription, or an elimination of a propagating mode.

## Smooth coefficient solution on the entire clock tube

The positive integrating factor mu=sqrt(-x)/w^(3/4) has mu_x/mu=A.
Consequently

    q(u,x)=mu(u,x)^(-1) integral[-1,x] mu(u,z)b(u,z) dz

is the unique solution with q(u,-1)=0. The integration segment stays
inside the open neighborhood stated above. Smoothness, and indeed
real analyticity locally, follow by differentiating a finite coefficient
integral on this regular domain. No assertion at x=0 is involved. The
independent audit uses the equivalent factor sqrt(-x)*p^(-3/2); their
ratio is the positive x-independent factor 2^(3/2)h^(3/4).

The bounds are uniform for every real u. The identity

    9h^2-h_u^2=9(1+u^2)^4(u^2-1)^2>=0

gives |h_u|<=3h. On the closed tube, |A|<=25/18 and |b|<=1/(4h^2).
Its maximal distance from the fixed base point is 1/10. The integral
Gronwall estimate and e^r<=1/(1-r), for 0<=r<1, give

    |q|<=e^(5/36)/(40h^2)<=9/(310h^2).

The parameter-dependent coefficient problem therefore has no finite-u
or tube-edge existence obstruction. At x=-1 its entire q(u,-1) function
vanishes; so do q_u(u,-1) and q_x(u,-1), since b(u,-1)=0. A value-only
clock check would nevertheless be insufficient to match the action on
the open tube, which is why the full ODE and scalar counterterm are kept.

## Exact scalar action and equations, with free matter retained

Set J2=F_repo(u,-x)+xq_u-3x(q+2f_u)^2/(16p^2). Then the literal
connection reduction has P=F_repo+xq_u, Q1=q, Q2=2q_x. At an arbitrary
normal frame the product rule gives

    partial_a x=2 phi^b phi_ab,
    div(q grad phi)=q Box(phi)+xq_u+2q_x phi^a phi_ab phi^b.

This tensor identity extends to any coordinates; the root algebra keeps
all ten H components, and the independent audit uses dense Fraction
fixtures in both signatures. It proves that the full lower action,
including F_repo away from X=1, is exactly the frozen target modulo
the displayed divergence. Dropping xq_u or treating c R_ab phi^a phi^b
as already in the scalar normal form would invalidate that conclusion.

The literal connection action is quadratic in the distortion and linear
in its derivatives before their boundary is removed. Projective gauge
fixing gives a nonsingular 60-variable algebraic stationary equation.
All 64 original Euler equations vanish on the chosen solution, including
the four dependent gauge equations. Thus the functional chain rule for
S[g,phi,chi,kappa_star(g,phi)] has no surviving connection-Euler term.
For compactly supported variations in the open domain, all boundary
terms are harmless. Its g,phi,chi equations are exactly the original
target's equations. Conversely a target solution lifts with kappa_star;
the full-rank theorem proves uniqueness modulo projective transformations.

This is stronger than a background reconstruction or principal-symbol
fit: it applies to arbitrary smooth target fields with timelike clock
gradient in the stated domain. The target's existing physical-g rolling
solution and free-chi backreaction are preserved. No new nonlinear
stability theorem or completeness statement for affine autoparallels
is inferred from that action equivalence.

## Quantitative algebraic map and units

In the trace-gauge component basis let E embed the 60 quotient variables
into the full 64-vector. The source projects as E^T j. Direct row sums
give ||E||_infinity=3 and ||E^T||_infinity=2. The independently factored
block inverses prove ||H_q^-1||_infinity<=880109/36000<25 throughout
the closed tube. Consequently the fixed-gauge stationary solution obeys

    ||kappa_star||_infinity <=150 ||j||_infinity.

These norms are in the orthonormal scalar rest frame and displayed
trace gauge; arbitrary Lorentz boosts are not claimed to preserve them.
The extra unsourced kernel at p^2=1/8 is outside the tube but is recorded,
not erased because the scalar forcing happens to be compatible there.

Under phi=tau*u, chi=M*chi_bar and physical coordinates=tau times their
dimensionless counterparts, X and the metric components are unchanged.
Curvature and two Hessians contribute tau^-2; one Hessian contributes
tau^-1. Hence p,c,J4 carry M^2, J3 and q carry M^2/tau, and J2 carries
M^2/tau^2. Four-volume gives overall M^2*tau^2. Delta=M^2 and the
physical quotient inverse norm bound is 25/M^2. This bookkeeping does
not introduce a kinetic term, heavy gap or physical cutoff.

The result is an exact classical auxiliary lift, not the adopted full
UV-matching gate. A kinetic connection completion would be a different
named action, with its own constraints, spectra, error estimates,
vacuum and finite-gravity analysis. Projective symmetry by itself is
not a proof of those properties.
