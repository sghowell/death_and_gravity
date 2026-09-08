# Exact ordered elimination, state and quantum boundaries

The complete source-centered vector action at fixed light fields is

    I_W=1/2 <W-S,M(W-S)>+zeta/2 <W,D W>.

The physical volume and index metric belong to the pairing. In the
timelike rest frame the algebraic mass bilinear has entries
(1/gamma_t,-1/gamma_s,-1/gamma_s,-1/gamma_s). D is the physical
Maxwell Euler operator, (D W)^mu=nabla_nu F(W)^(nu mu). Integration
by parts gives <W,D W>/2=-integral F(W)^2/4, plus its specified
surface term. Neither M nor D is replaced by a scalar constant; M
does not generally commute with D. The affine complement and complete
source come from the all64 frozen elimination.

Let K=M+zeta D. On a common domain with an appropriate inverse,

    W_star=K^-1 M S=S-zeta K^-1 D S,
    M-M K^-1 M=zeta D-zeta^2 D K^-1 D.

The latter follows by substituting M=K-zeta D, with no interchange
of factors. Symbolic noncommuting operator identities and independent
dense noncommuting symmetric matrices check the ordering.

For a compatible symmetric variational boundary prescription, complete
the square exactly:

    I_W=1/2 <W-W_star,K(W-W_star)>+I_induced,
    I_induced=zeta/2 <S,D S>-zeta^2/2 <D S,K^-1 D S>.

No power series or convergence hypothesis was used. The first local
term is the full source curl density bounded in source.md. Its exact
remainder still contains the differential inverse. In a specified
positive norm, an independently established bound ||K^-1||<=L would
imply a corresponding bilinear remainder bound zeta^2 L||D S||^2/2,
with the pairing controlled in that norm. This checkpoint supplies
no such general Lorentzian operator estimate and makes no numerical
claim about that nonlocal remainder. Neither a rest-frame algebraic
mass bound nor a background canonical-frequency floor supplies it.

## Retarded dynamics is not a one-copy stationary action

For a compatible initial-value Green operator, the full linear vector
equation at fixed light fields can instead be written

    W=S-zeta G_ret D S+W_hom.

The homogeneous part must be fixed by the physical Cauchy data and
temporal constraint. It may be removed only after their compatibility
with the source has been checked. This response formula is not an
instruction to replace K^-1 by G_ret in the stationary-action formula.

Varying 1/2<J,G_ret J> in a single-copy action gives
(G_ret+G_ret^T)J/2, not G_ret J. A two-time lower-triangular exact
control has no late source in its early retarded row, but does have
one in the action's early variational derivative. This directly
detects the advanced contribution. A doubled initial-value/in-in
formalism, or the full original constrained equations, is needed
for causal feedback. S6.43's response theorem remains a state-specific
result on the rolling coefficients, not a variational action theorem.

The same issue is illustrated in
[Galley, arXiv:1210.2745](https://arxiv.org/abs/1210.2745), equations
(2)-(4). The present ordering and finite-dimensional causal controls
are derived directly; no solution or bound for this model is imported.

## The source alignment does not remove quantum loops

At a finite regulator and a stated Gaussian contour/measure, completing
the square leaves a vector determinant det(K)^(-1/2). Its formal
in-out effective-action contribution is (i hbar/2) Tr log K; a Euclidean
form requires its own valid continuation. K is independent of S but
depends on the light metric and clock through M and D. Thus S=0 or
vanishing first source variation does not make the vector determinant
constant or zero. The all64 auxiliary measure, field transformations,
light loops, subtractions and finite counterterms are separate issues.

[Ruf and Steinwachs, arXiv:1806.00485](https://arxiv.org/abs/1806.00485),
sections V-VI, distinguish the nonzero massive Proca determinant from
the gauge-field case and treat a general background mass tensor.
Their analysis supplies a relevant route for the next quantum
calculation, not an already evaluated determinant or error bound
for the current rolling parent. No quantum, cutoff or UV verdict
follows from the classical Gaussian identity. Original P8 stays open.
