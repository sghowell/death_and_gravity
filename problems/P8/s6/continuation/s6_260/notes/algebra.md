# Full Grassmann jet algebra and canonical boundary

Use a left odd derivation s and the graded Leibniz rule
s(AB)=(sA)B+(-1)^parity(A) A(sB). Spatial c^i and all their derivatives
are anticommuting. The clock is fixed, c^0=0, but partial_0 c^i is not zero.
Set sc^i=c^j partial_j c^i.

For any spatial tensor density T of weight w, sT=Lie_c T. Every lower
index contributes +(partial_i c^j)T_...j..., every upper index contributes
-(partial_j c^i)T^...j..., and the density contributes w(div c)T.
The explicit scalar, covector, vector-density, symmetric metric,
metric-momentum and Q rules are retained. Q has weight 2/3; canonical
momenta have weight one. Scalar rules apply separately to M1, H, N
and T, without altering any coefficient function of the fixed clock.

These rules satisfy s^2=0 because sc supplies the ghost self-commutator,
the second-derivative terms cancel by graded antisymmetry and ordinary
jet derivatives commute. The code expands every exterior coefficient
on independent exact bosonic jets. Symmetric metric/momentum restrictions
have six components each. The whole affine connection has 64 components,
not a torsionless subset. For its derivative-index-last convention,

s Gamma^a_bc =
 c^j partial_j Gamma^a_bc
 -Gamma^d_bc partial_d c^a
 +Gamma^a_dc partial_b c^d
 +Gamma^a_bd partial_c c^d
 +partial_b partial_c c^a.

Here d runs over all FOUR spacetime indices. Even though c^0=0,
the d=0 upper-index term with partial_0 c^a generally survives.
Independent variation of the full flat-metric Levi-Civita connection
reproduces every inhomogeneous second-ghost-derivative component.

For the shift multiplier,
sN^i=partial_0 c^i+c^j partial_j N^i-N^j partial_j c^i.
The full h=C gamma ADM map then transforms as a four-metric. The
original W_0=N T+N^i W_i transforms as a four-covector, including its
partial_0 c^i W_i contact.

## Canonical action and regular branch

Use the entire S258 spatial generator H_i, including metric-density
terms, vector curl AND -W_i div pi_W, both matter pairs and the N,T
primary pairs before restriction. For the canonical density

L=pi^ij dot gamma_ij+pi_W^i dot W_i+sum p_A dot f_A-H0-N^i H_i,

H0 is the entire S257 scalar-density Hamiltonian, not a replaced
potential. Its scalar-density property follows from that full action.
Direct graded differentiation gives

sL=partial_j(c^j L+F^j),
F^j=2 pi^jk gamma_ki dot c^i+pi_W^j W_i dot c^i.

The extra flux is not zero pointwise. Periodicity or the stated spatial
boundary conditions remove only its integral. The original primitive
boundary and endpoint phases are retained. Restricting the two
auxiliary primary pairs leaves eleven retained coordinates and
twenty-two phases before spatial gauge reduction.

The unique local auxiliary branch is equivariant by the S257/S258
regularity and scalar-density argument, so the geometric rules induce
the classical reduced BRST differential there. They are not asserted
to equal an off-constraint Dirac/BFV differential on every eliminated
lapse, vector-normal, jet or shift-primary variable.

## Complete off-gauge fermion

Let chi^i=partial_j Q^ij and s bar c_i=B_i, sB_i=0. With fixed gauge
width alpha, Psi=integral bar c_i(chi^i+alpha B_i/2) gives

sPsi=integral [B_i chi^i+alpha B_i B_i/2-bar c_i M^i_j c^j].

M is exactly the full S258 operator, including all chi and derivative
of chi terms BEFORE restriction. Direct exterior algebra verifies
s(Mc)=0 and s(sPsi)=0. This is a classical gauge action, not a proof
that a separately regulated quantum measure is invariant.
