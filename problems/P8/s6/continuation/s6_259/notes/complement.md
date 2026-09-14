# Complete 56-direction determinant and exact inertia

Start from the actual 60-dimensional quotient Hessian M of the original
unrestricted affine connection. Let N be the full four-trace map and
D=N M^-1 N^T. The source-pinned S174 update is

    Mnew=M+N^T Xi N,        Xi=eta-D^-1,
    L=M^-1 N^T D^-1,       N L=I,

with eta=diag(1,-1,-1,-1). No four-vector ansatz is used to obtain M.

Choose the pivot columns (0,1,4,8) of N, make their right inverse E,
and let the remaining 56 free-coordinate columns be U. Set

    K=U-E N U.

Then NK=0, the free-coordinate submatrix of K is I, and [K,E] has
determinant -1. Since N(L-E)=0, L-E=K C for some C, so replacing E
by L is a unit block shear. Therefore det[K,L]=-1 everywhere.

The complete transformed quadratic form is

    [K,L]^T Mnew [K,L] = diag(A,eta),        A=K^T M K.

All 56-by-4 mixed entries and all 56-by-56 mass-update differences
are checked. The inverse is assembled from every actual block of A,
not inferred from a trace determinant alone.

The nine block sizes are 6,6,8,6,8,8,8,3,3. Their determinants, in
that order, are

    p^4(8p^2-1), p^4(8p^2-1), -96p^5(2p^3-1),
    p^4(8p^2-1),
    1024p^12(8p+5), 1024p^12(8p+5), 1024p^12(8p+5),
    -16p^6, -16p^6.

Their product gives the determinant in the formulation. Independently,

    det A = det[K,L]^2 det M det D,
    det M = -4503599627370496 p^72 (8p^2-1)^3.

This Schur identity agrees with the separately computed nine-block
product. It is not an omitted determinant factor.

For p>0 and 1/8<p^2<3/10, 8p^2-1>0 and 8p+5>0.
Also p<3/5 implies 2p^3<54/125<1. Thus det A>0 and no eigenvalue
is zero in the connected strict interval. Neither endpoint is included.

At p=1/2 the report gives complete rational congruence transformations
for all nine blocks. The algorithm uses nonzero one-dimensional
diagonal pivots or a two-dimensional off-diagonal pivot with negative
determinant; the latter contributes one positive and one negative
direction. Every congruence is checked by multiplication. The block
inertias sum to (26,30), not (56,0).

A continuous real symmetric nonsingular matrix cannot change inertia
without an eigenvalue passing through zero. Hence (26,30) holds on
the whole strict interval. Floating eigenvalues at several rational
points are independent diagnostics only; they are not the proof of
this all-parameter assertion.

The indefinite directions are algebraic auxiliaries, not by themselves
propagating ghosts. Their role in a finite oscillatory integral and its
canonical constrained density is a separate calculation.
