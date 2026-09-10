# Dirac norm and the full regulated Ward sum

Use Hermitian Euclidean gamma matrices with
{gamma_mu,gamma_nu}=2 delta_mu_nu. The explicit four-by-four
representation is checked entrywise. For real Euclidean q,

    S0(q)=(m-i gamma.q)/(q^2+m^2),
    S0.dagger S0=I/(q^2+m^2).

Thus ||S0||=(q^2+m^2)^(-1/2). A real unit transverse
polarization gives ||gamma.epsilon||=1. A complex shift
obeys ||gamma.P||<=sum_mu |P_mu| by the triangle inequality,
because each gamma matrix has operator norm one.

The exact resolvent identity is

    S(q) i gamma.k S(q+k)=S(q)-S(q+k).

Its cleared numerator is checked for independent components;
as a polynomial identity it also holds at complex k.
Insert a contracted gauge vertex in the three positions of
each of the two cyclic orders of the remaining vertices.
The adjacent triangle differences telescope, using
translation-invariant regularization. The scalar vertices
commute with the color generator. The remaining nonabelian
commutator has trace zero, as independently checked with all
64 SU(3) generator pairs. A possible triangle connected to a
three-gauge vertex also vanishes by its single-generator trace.

Apply these statements to the FULL six-box sum before
removing the regulator. Each individual zero-momentum box
can contain a logarithmic UV divergence; it is not assigned
a finite four-dimensional norm integral. The Ward identities
forbid a neutral-scalar Phi^2 A^2 operator, so the complete
degree-zero tensor vanishes. There is no degree-one
parity-even rank-two Lorentz tensor: it would require an
invariant rank-three tensor. Real scalar Yukawas do not
supply a gamma5 or CP-odd coefficient.

For degree two, use three independent external vectors
k,l,p. The complete parity-even tensor basis has eta times
their six Gram products and their nine ordered outer
products. Contract with k on the first gauge leg and l
on the second, and equate every independent Gram coefficient.
The resulting 36-by-15 rational matrix has rank 14.
Its one-dimensional nullspace is

    (k.l) eta_mu_nu - l_mu k_nu.

This classification is off shell, before restricting the
light mass. Derivatives on the neutral scalar cannot add an
independent gauge-invariant operator at this degree.
The finite coefficient is therefore fixed by the constant-
scalar heavy determinant already derived in S6.128:

    C Phi^2 F^a_mu_nu F^a_mu_nu,
    C=-a Y/(48 pi^2 m^2)=-a Y/(3 Q m^2),
    a=g_s^2, Y=y^2, Q=16 pi^2.

The degree-two terms and all higher degrees are individually
UV convergent. Their regulator limit equals the ordinary
four-dimensional integral; no evanescent factor times an
unremoved degree-zero pole is silently discarded.
The regulator input and its nonanomalous restriction are
identified in literature.md. This is a written field-theory
argument with exact matrix and tensor checks, not a formal
proof of the regulator construction.
