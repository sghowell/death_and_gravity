# Ordered adiabatic marker and physical-current parity

The formal marker zeta is distinct from the physical amplitude epsilon.
For the same finite tenth-order reference, write

    rhat(zeta)=sum_(n=1)^10 zeta^n r_n.

The exact ordered recurrence is
r1=-S/(2i omega) and
r_(n+1)=(r_n'-[R,r_n]+sum_(j+l=n)r_j S r_l)/(2i omega).
Real R,S,omega and real time derivatives imply inductively that r_n
is i^n times a real symmetric matrix. The transpose of the commutator
has the correct sign because R is skew, and the sum over both ordered
pairs j,l preserves symmetry. No simultaneous diagonalization is used.

For complex zeta define rsharp(zeta)=sum zeta^n r_n^dagger; do not
conjugate zeta. Put C=(I-rsharp r)^-1 and use the coefficient adjoint
Fsharp of F=[I+r;-i(I-r)]/sqrt(2). Extend the physical real covariance
holomorphically by averaging W=F C Fsharp with its coefficientwise
complex conjugate. At real zeta this equals Re W and hence the actual
finite-reference covariance.

Every total-order n coefficient in the QQ and PP blocks is i^n
times a real matrix; their odd coefficients vanish after averaging.
The QP/PQ blocks have an additional factor i, so their odd terms
generally do not vanish. The physical G_D is real block diagonal,
therefore tr(G_D Sigma(zeta)) is even. The same argument holds for
the first two physical amplitude derivatives since they preserve the
real-phase property and the block structure.

The code reconstructs the full noncommuting reference through order10,
checks its phase and symmetry at every order, and checks complete
covariance parity through order6. A negative control retains a
nonzero first-order QP block: discarding the whole odd covariance
would be wrong even though the metric current is even.
