# The original trace gauge needs spatial compensation

Keep F_mu=Gamma^a_{a mu}, the original four projective columns G, and
the source-pinned relation F G=4I. The independent projective parameter
acts as delta Gamma^a_bc=delta^a_b eta_c.

Taking the full trace of the spatial connection rule gives
sF_mu=Lie_c F_mu+partial_mu(div c)+4 eta_mu.
The inhomogeneous term is present even for F=0. Therefore spatial
diffeomorphisms alone do not preserve this gauge. On the slice use

eta_mu=-partial_mu(div c)/4,

including its time component. The denominator is four, not the
number of spatial dimensions. The compensated 64-component
connection rule remains nilpotent.

Before fixing the projective gauge introduce all four odd eta_mu,
with s eta_mu=Lie_c eta_mu, interpreting both factors with the
graded left rule. The full spatial/projective semidirect differential
is nilpotent. Set g_mu=F_mu/4 and

zeta_mu=sg_mu=eta_mu+Lie_c g_mu+partial_mu(div c)/4.

Then szeta_mu=0. Thus g,zeta are contractible pairs. The change
eta to zeta is triangular with unit finite Berezin Jacobian; it is
not permission to ignore state, source or continuum regulator terms.

## Whole mixed ghost matrix

Before restriction the seven-generator gauge derivative has blocks

[ 4 I_4,  Lie_xi F + d div xi ]
[   0,                M(Q)    ].

The finite block determinant is 256 det M. Normalizing the trace
gauge to F/4 changes the projective block to the identity. Neither
convention changes the physical quotient. The proof is a finite
triangular matrix identity, not multiplication of separately
regularized continuum determinants.

At flat Q=I, nonzero spatial k and time frequency omega,
M=k^2 I+kk^T/3 and the mixed block is
-(omega,k_1,k_2,k_3)^T k^T. The complete 7-by-7 matrix and both
inverse products are checked, including the omega row. The
determinant is (1024/3)(k^2)^3, independent of omega.
No inverse on the residual translations is claimed.

## Original source is retained

The complete original source obeys G^T J=0. The exact dual
projector (I-GF/4)^T constructs such a 64-component source, and
its pairing with the compensating projective shift is zero.
This removes only that gauge shift. It does not remove the
nonlinear composite connection map, the S259 inverse quadratic
source contact, or their derivatives and physical boundary terms.
