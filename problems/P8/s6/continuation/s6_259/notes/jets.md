# Finite derivative-coordinate lift and constrained density

The actual source-centered connection shift contains metric and clock
derivative data. It is therefore not an ordinary point transformation
on the unextended field configuration space.

## Retain the defining derivative constraints

For a finite derivative order, introduce independent coordinates for
those derivatives. For example introduce v for qdot and a multiplier
p imposing qdot-v=0. A higher derivative uses a finite chain
v_next-d(v_previous)/dt=0. Spatial derivatives at a common finite
regulator can likewise be retained as their exact prescribed
constraints; this does not assert a diffeomorphism-preserving regulator.

The enlarged first-order action is the original action evaluated at
the independent derivative coordinates plus all these kinematic
multiplier terms. Integrating those multipliers enforces the original
relations. No additional physical derivative mode has been introduced.
Endpoint data and integration-by-parts boundaries remain prescribed.

On this extended space the whole invertible connection change is a
point map x=F(Q). Its cotangent lift is

    P=DF(Q)^T p_old,        p_old=DF(Q)^-T P.

It preserves the full canonical one-form and symplectic form. Its
phase-space Jacobian is one even if its configuration determinant
is not. Every derivative-coordinate momentum shift is included.
The report independently checks an eight-phase nonlinear
velocity-dependent example. The general chain rule proves the result
for all coordinates of the actual 64-component map, not just those
of that diagnostic.

## Eliminate only the regular algebraic pairs

After the exact source-centered change, the complement part of the
action is lambda z^T A(Q) z/2. Here Q includes the independent
derivative coordinates and the retained fields; lambda and A may
depend on them. Its auxiliary primary momenta are pi_z=0.
The 56 secondary equations at zero complement source are lambda A z=0,
up to the common Hamiltonian sign convention. Since A is nonsingular
there is one root z=0.

For a general regular source-dependent version write the secondaries
as C(Q,P,z)=0 and D=partial_z C. The full auxiliary bracket is

    M_z = [[0,-D],[D^T,E]],

with E the full secondary bracket, not silently set to zero.
Its inverse has lower-right block zero:

    M_z^-1 = [[D^-T E D^-1,D^-T],[-D^-1,0]].

Thus det M_z=(det D)^2. The positive second-class density
sqrt(det M_z)=|det D| cancels the delta(C) integration Jacobian
on the single regular root. Retained functions independent of
z,pi_z commute with the primaries, so their Dirac bracket is
the original canonical bracket. This is a finite-dimensional
identity, applied to every cell of a common finite regulator.

The remaining kinematic, metric, lapse and spatial constraints are
not discarded. They are pulled back to the auxiliary root and use
the induced Dirac bracket. One can see the same statement for a
further finite constraint matrix by block Schur factorization:
its full determinant is det M_z times the determinant of the
remaining Dirac-bracket matrix. This does not require the remaining
matrix to be invertible before gauge treatment, and it does not
by itself solve the residual constraints or specify a quantum state.

The projective coordinates have their four primary momenta and the
specified algebraic gauge. They are not dynamical vector modes.
The four retained trace/Proca directions remain in the physical
retained system with their complete nonlinear source.

## An invalid shortcut and its exact counterexample

Use Lred=(v^2-q^2)/2, f=(qv+q^2,v^2-q) and

    A=[[2+q^2,qv],[qv,-3-v^2]],
    L=Lred+(y-f)^T A(y-f)/2.

A is nonsingular and indefinite at every real q,v. Eliminating y
classically gives the regular free oscillator. Yet the unextended
velocity Hessian at y=f is

    L_vv|y=f = 1+f_v^T A f_v,

and at q=0 it vanishes when v^2=(sqrt(10)-3)/2.
Therefore an argument that globally inverts that unextended Hessian
would fail. The finite-jet formulation instead retains v and its
kinematic constraint, reduces the regular algebraic pairs, and
then gives exactly Hred=(p^2+q^2)/2. All six primary/secondary
brackets in this diagnostic are retained and checked.

## Quantum ordering is an additional prescription

These are classical cotangent and finite constraint-density facts.
A finite conditional Gaussian insertion can lift an already
specified retained functional without altering its weight.
Neither assertion licenses replacing the transported time-sliced
action, operator prescription or initial state by a new one.

For example, the free particle under x=exp(Q), with wavefunctions
transformed as half-densities, has

    p_transformed=-i hbar [exp(-Q) partial_Q-exp(-Q)/2].

Its exact squared Hamiltonian differs from the naive Weyl
quantization of P^2 exp(-2Q)/2 by
hbar^2 exp(-2Q)/8. The report derives the full differential-operator
identity. This nonzero term prevents inferring complete nonlinear
quantum equivalence from a unit Liouville Jacobian alone.
