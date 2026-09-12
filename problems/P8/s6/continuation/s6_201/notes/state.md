# Initial occupation weights, without state momentum derivatives

The actual all-order prepared state has the source-pinned bound

    |beta0(k)| <= (B6+E/m^4) nu^-6 < B nu^-6,
    B6=1813229, E=5347035781757616, B=2e6, m=1000.

This follows already from the uniform mixing envelope at the
preparation surface. No derivative of an initial cutoff function or
Borel coefficient is taken. In particular |beta0|^2<=B^2 nu^-12<1.

The normalized exact state satisfies |alpha0|^2-|beta0|^2=1.
For every full creation pair, alpha_k alpha_l multiplies both detector
and source amplitudes. Their product with one conjugated amplitude
therefore contains |alpha_k|^2|alpha_l|^2, independent of both initial
phases. Each contact covariance contains |alpha_k|^2.

Writing b_k=|beta_k|^2 gives the exact increment
b_k+b_l+b_k b_l. The quadratic term is retained. Since 0<=b_k,b_l<1,

    2(b_k+b_l)-(b_k+b_l+b_k b_l)
      =b_k(1-b_l)+b_l>=0.

This proves the stated complete occupation majorant. The actual
phases and mixing are not reset; only their effect on the algebraic
comparison current is evaluated.

The unit-W8 canonical pair f=(2W)^(-1/2) exp(-i integral W),
p=f'-d f has f conjugate(p)-p conjugate(f)=i for real W>0 and d.
Multiplication by alpha0 changes this to i(1+b_k). Removing that
prefactor restores the Wronskian, not the exact field equation:
W8 still has its original nonzero residual. A direct independent
variable-W oscillator counterexample checks that distinction.
