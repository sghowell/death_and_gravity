# Canonical matter boundary and full seven-mode quadratic replay

The canonical variables used in the nonlinear Hamiltonian
must be distinguished from those of the earlier quadratic
scalar action. A mismatch here changes the apparent
constraint source and cross terms.

## Both matter and curvature momenta shift

Let P0=a^3 ell=1/10 in global clock coordinates. Expanding
the original canonical matter term gives the quadratic
cross term 3P0*v*chi_dot. The earlier scalar action writes
it instead as -3P0*chi*v_dot, with their difference

d_t(3P0*v*chi).

The exact background identity P0_dot=0 matters. Thus the
old quadratic momenta are related to the current spatial
canonical variables by

p_old=p-3P0*chi,
Pchi_old=Pchi-3P0*v.

Both shifts are required. Their one-form difference is
the differential of 3P0*v*chi; the two canonical pairs
and their vanishing cross momentum bracket are checked.
No additional explicit time term remains from this
boundary because P0 is constant. This comparison does
not redefine the nonlinear fields of the old action.

The leading York source p-3P0*chi is therefore exactly
p_old. In local point units a=1 this is the matter shift
found independently from the actual momentum generator.

## Full spatial geometry and actual lapse correction

For each labelled pair compute the inverse metric,
volume, Christoffels and scalar curvature from the full
linear spatial gauge, not just its conformal part.
Solve the full momentum constraints to the required
degree and form all nine S6.75 invariants with their
actual geometric contractions.

In particular, if piMixed=pi*g/sqrt(g), the trace
invariant is dp=(2/3)tr(piMixed)+2H and the matter
invariant is dc=(ell+Pchi)/sqrt(g)-ell. At first order

dp=p/3, dc=Pchi-3ell*v, R3=4q*v.

The full Hamiltonian gives, by direct lapse derivatives,

Hcal_N,dp=3Theta,
Hcal_N,dc=-w,
Hcal_N,j=0,
Hcal_N,R3=-Lambda/2.

These are independently computed from its actual
coefficient functions. The source alignment is responsible
for the zero linear vector-divergence coefficient.
The boundary primitive and margin remain included in
the background Hessian Hcal_NN=-2J_e.

It follows that the first-order lapse force is

F1=Theta*p-w*Pchi+(3ell*w-2Lambda*q)*v.

Substitution of the actual stationary lapse contributes
F1^2/(4J_e) to the quadratic normal Hamiltonian.
This uses stationarity of the full action; it does not
freeze its coefficient derivatives or introduce an
inverse Theta.

## Moving canonical terms are kept

The exact linear spatial gauge uses the gravitational
one-form boundary of the original general York construction.
At a local evaluation point its Hamiltonian contributions
are

-2H*pi:g
-(H_dot+3H^2)*(6v+3v^2-gamma:gamma/2).

The shift from the actual background matter field adds
-ell*Pchi, apart from an irrelevant field-independent
background term. P0 is the canonical momentum density;
replacing it by a time-dependent unweighted ell before
this transformation would produce a wrong boundary.

The full volume multiplying the invariant Hamiltonian,
the spatial curvature, all these moving terms and the
lapse correction are kept before extracting degree two.
In particular H_dot is not set to zero at the bounce.

## Compact scalar result and vector/tensor blocks

Write q=|k_physical|^2 in the fixed local bounce-time
units, and use the current canonical momenta p,Pchi.
The scalar/matter Hamiltonian is

H2=(Pchi-3ell*v)^2/2-ell*p*chi/2
   +(q/2+3ell^2/4)*chi^2-q*v^2+R^2/J_e,

R=-Theta*p/2+(Lambda*q-3ell*w/2)*v+w*Pchi/2.

Direct substitution of the two momentum shifts into the
original S6.56 regular scalar Hamiltonian gives this same
expression. All scalar phase pairs are checked against
the full spatial Hamiltonian, not merely against this
boundary identity.

The tensor block is
2*Pi_TT:Pi_TT+(partial gamma_TT)^2/8.
Both independent TT polarizations are retained with
their dual canonical momentum tensors.

The vector block is the actual Proca Hamiltonian

Pi^2/(2zeta)+(partial_i Pi^i)^2/2
+zeta*F_ij F_ij/4+W_i W_i/2.

It has all three spatial polarizations. There is no
quadratic scalar/vector or tensor/vector source, even
though the vector participates in nonlinear momentum
constraints and higher interactions.

The implementation uses fourteen phase channels:
two scalar canonical pairs, two tensor canonical pairs
and three vector canonical pairs. Every one of the
105 independent phase-pair coefficients is extracted
at exact symbolic u and positive wave scale and compared
with the prior scalar/matter and tensor/vector blocks.
For a real integrated quadratic Hamiltonian the remaining
entries obey H_ji(k)=H_ij(-k); the reported full matrix
uses that functional-Hessian identity, not an assumed
point-frozen self-adjoint coefficient matrix.

All these comparisons are polynomial/rational identities
on the declared nonzero-momentum domain. They do not
diagonalize a resummed quantum symbol, drop growing
residues, or derive a new physical cutoff.

## Relation to higher interactions

The canonical coordinates in this checkpoint are the
ones to use for the actual cubic/quartic reduction.
The quadratic match checks that the geometric, source
and moving-boundary conventions are consistent.
It does not itself supply those higher vertices.
S6.75's invariant Taylor bounds must still be composed
with this constrained geometric map, with the declared
hard-transfer and finite-time observable domain.
