# Whole conformal density gauge and its nonsymmetric symbol

For positive gamma define Q^ij = det(gamma)^(1/3) gamma^ij.
It is a contravariant tensor density of weight w=2/3, not an ordinary
tensor and not the weight-one harmonic density.

For the complete physical map h_phys=C(u,N) gamma with C>0,

    det(h_phys)^(1/3) h_phys^-1
      = (C^3 det gamma)^(1/3) C^-1 gamma^-1 = Q.

This pointwise identity holds for spatially varying N and C. Taking
spatial derivatives preserves it; no lapse-gradient terms have been
dropped. With gamma=a_hat^2 exp(2v) exp(t), trace t=0, Q=exp(-t).
The nonlinear spatial gauge is chi^i=partial_j Q^ij=0.

The exact infinitesimal density variation is

    delta Q^ij = xi^k partial_k Q^ij
                 - Q^ik partial_k xi^j - Q^jk partial_k xi^i
                 + w Q^ij partial_k xi^k.

Apply partial_j before setting chi to zero. Product-rule terms with
one derivative of Q and one derivative of xi cancel in pairs after
renaming the contracted indices. The result is

    M xi = xi.grad chi - chi.grad xi + w chi div xi
           - Q^jk partial_j partial_k xi^i
           + (w-1) Q^ij partial_j div xi.

Thus the gauge-surface-only operator is not the full off-gauge
derivative. Its omitted chi terms can contribute to nonlinear vertices.
They are kept throughout the Fourier construction.

For a real covector k != 0 the whole principal symbol is

    sigma_M(k) = r I + (1/3)(Qk) k^T,     r=k^TQk.

Every vector in ker k^T has eigenvalue r; Qk has eigenvalue 4r/3.
The determinant is 4r^3/3. Although the symbol is generally not symmetric
in the Euclidean metric,

    Q^-1 sigma_M = r Q^-1 + kk^T/3

is symmetric positive definite. This supplies the principal symmetrizer.
It does not by itself prove coercivity for an arbitrary variable
coefficient field; the full weak-form estimate is in the local note.

At Q=I the nonzero-mode inverse is

    M0(k)^-1 = [I-kk^T/(4k^2)]/k^2.

The inverse is not assigned at k=0. Constant vectors lie in the kernel
on the exact slice, whereas off the slice their action is xi.grad chi.
Keep that distinction in nonlinear variations and residual-group scope.
