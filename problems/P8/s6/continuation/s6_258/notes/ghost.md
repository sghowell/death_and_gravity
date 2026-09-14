# Entire off-gauge Fourier operator and finite ghost vertices

Let Q_ell be a symmetric shape-density Fourier coefficient, k the
incoming ghost momentum and p=k+ell the outgoing momentum.
Fourier transforming the complete density Lie derivative gives

    M(p,k) = (p^T Q_ell k)I + (Q_ell k)p^T
             - (Q_ell p)ell^T - (2/3)(Q_ell p)k^T.

All three ghost components remain. The outgoing p=0 row vanishes
identically because chi is a divergence. The incoming k=0 column
is -(Q_ell ell)ell^T. It is not generally zero off the exact gauge
surface. When Q_ell ell=0, the full expression reduces to the
gauge-surface symbol (k^T Q_ell k)I+(Q_ell k)k^T/3.
No such simplification is made before off-gauge differentiation.

## A complete exponential shape family

Take S=diag(0,1,-1) and gamma=a_hat^2 exp(2v) exp(t S cos x1).
Then

    Q=diag(1,exp(-t cos x1),exp(t cos x1)),

with det Q=1 and div Q=0 exactly at every t, not just through some
finite order. For |t|<=1/100,

    ||Q-I||op <= exp(|t|)-1 <= |t|/(1-|t|) <= 1/99 < 1/8.

The local weak-coercivity estimate therefore applies. Project the
entire spatial operator to the eight modes (a,b,0), a,b in {-1,0,1},
excluding (0,0,0), retaining all three ghosts: a 24-by-24 matrix.
The reference and first two exact derivatives have coefficient data

    Q0(0)=I,
    Q'(e1)=Q'(-e1)=-S/2,
    Q''(0)=S^2/2,
    Q''(2e1)=Q''(-2e1)=S^2/4.

Write A=M(0), B=M'(0), C=M''(0). The report retains the entire matrices
and det A, not only a selected polarization. Differentiating
M M^-1=I gives (M^-1)'=-A^-1 B A^-1. At finite dimension Jacobi's
identity gives

    (log det M)' = Tr(M^-1 M'),
    (log det M)'' = Tr(M^-1 M'' - M^-1 M' M^-1 M').

The relative determinant det M(t)/det A retains the explicitly
specified reference normalization. Its first derivative is zero and
its second logarithmic derivative is exactly 3. The separate value
obtained by deleting C is displayed and differs from 3. Therefore
keeping only two first-vertex insertions is wrong even in this
finite gauge-satisfying example.

Independent numerical tests use the entire Fourier coefficients of the
exponentials, expressed with modified Bessel functions I_n(t), before
projecting. The yy coefficient is (-1)^n I_n(t), the zz coefficient
is I_n(t), and xx is present only at n=0. Central differences of the
full matrices and their relative determinants recover B,C and 3.
Three-dimensional nonzero-mode subspaces independently test the full
nonsymmetric weak-coercivity margin.

## Scope

The finite determinant is a spatial diagnostic, not a spacetime
regularized path integral. A canonical Faddeev-Popov prescription would
also need the momentum constraints, gauge delta functions, residual
group treatment, contours, time regulation, state and a regulator with
the appropriate quantum identities. None is produced by merely
normalizing this finite determinant. No continuum ghost subtraction,
vanishing of quantum determinants, BRST regulator, covariance of the
full quantum measure or physical loop size is inferred.
