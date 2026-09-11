# Complete physical metric-vertex transfer

Keep the S6.188 Hamiltonian M=diag(V,K), Z=(A,pi), and the S6.190
balanced frame Z=T Y, T=diag(B/sqrt(omega),sqrt(omega) B^-1).
Here B is the actual positive matrix square root of K, including its
longitudinal constraint term. KV=omega^2 I implies T^t M T=omega I.

Let D(t) be any fixed smooth real symmetric tracefree matrix detector
with operator norm at most one. It differentiates the shear coordinate
gamma, not a scalar proxy. The actual conjugate metric current is

    J_D=-tr(M_D C_phys)/2=-omega tr(G_D Sigma)/2,
    G_D=T^t M_D T/omega.

All bounds are uniform in direction of momentum and in this detector.
The volume factor a^3 is independent of the unimodular shear.

For any n<=3, the full matrix-exponential Frechet derivative is the sum
over all n! ordered simplex integrals. Its norm is at most
exp(delta) times the product of the direction norms: the simplex
volume is 1/n!, and the exponential segment lengths sum to one.
No matrix products are commuted. Comparing with exp(+/-gamma)>=
exp(-delta)I gives, separately for the positive mass, electric and
magnetic pieces, a quadratic-form bound exp(2delta) times the base
piece. Adding the pieces gives

    ||M^-1/2 M_(D Gamma^a) M^-1/2|| <= exp(2delta)<2, a=0,1,2.

The constraint term has zero shear derivatives but remains in M.
This argument does not require its eigenvectors to commute with D.

Write p_j=omega_(epsilon^j)/omega. S6.191 gives |p1|<=1, |p2|<=2,
||B^-1 B_e||<=4, ||B^-1 B_ee||<=68,
||(B^-1)_e B||<=4 and ||(B^-1)_ee B||<=100.
Transposing the last two products gives the bounds in the momentum
block. Full raw product differentiation therefore yields

    ||T^-1 T_e|| <= 4.5 < 5,
    ||T^-1 T_ee|| <= max(68+4+1.75,100+4+1.25) < 110.

These are bounds on raw second derivatives, not on the derivative of
the logarithmic first-derivative matrix.

Let H_j=T^-1 T_(epsilon^j), and let G_(D Gamma^a) mean the normalized
physical higher Hamiltonian vertex, without differentiating T.
For F=T^t M_D T, all ordered product terms give

    F_e/omega = H1^t G + G H1 + G_(D Gamma),
    F_ee/omega = H2^t G + G H2 + 2 H1^t G H1
                 + 2 H1^t G_(D Gamma) + 2 G_(D Gamma) H1
                 + G_(D Gamma Gamma).

Their bounds are 22 and 582. Differentiating 1/omega as well gives

    ||G||<2, ||G_e||<24<25, ||G_ee||<634<650.

The code verifies both complete formulas against independent raw
polynomial products with noncommuting matrices. The proof supplies
the continuum operator bounds; fixtures do not replace it.
