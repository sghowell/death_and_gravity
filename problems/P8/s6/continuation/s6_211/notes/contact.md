# Full contact cancellation only in the spatial difference

The exact unimodular Hamiltonian is

    H=1/2 integral[pi^T exp(gamma)pi/a
       +(div pi)^2/(a^3m^2)
       +a m^2 A^T exp(-gamma)A
       +(curl A)^T exp(gamma)(curl A)/a].

The temporal constraint A0=-div pi/(a^3m^2) is gamma independent. For independent noncommuting tracefree perturbations D,Gamma the mixed exponential coefficient is

    H_DGamma=(D Gamma+Gamma D)/2.

Consequently the original complete second vertex is

    M_DGamma(k,q)
      =diag(a m^2 H_DGamma+Ck^T H_DGamma Cq/a,H_DGamma/a).

The positive mass contact is essential; the first-variation negative mass sign cannot be reused here.

For a homogeneous covariance, the quadratic-current trace is diagonal in internal momentum. The detector/source Fourier factors pair as D_hat(-P),Gamma_hat(P); their pointwise mixed convolution has total0. Therefore q=k in the second vertex, with no remaining transfer P. This is exact before integrating internal k or applying the original one-leg cutoff chi_K(k). It holds for the actual covariance and the unit-W8 reference covariance separately.

Thus C_K(P)-C_K(0)=0. This is not C_K(P)=0: for a simple diagonal covariance at k0 the current is

    -(a m^2 c_A+c_pi/a) tr(D Gamma)/2,

which is generally nonzero. Noncommuting tensor tests check both the exact cancellation and a nonzero full trace. The magnetic, electric, mass and constraint contents are all retained.

For smooth compact spacetime smears the same argument applies inside the Fourier pairing. P0 denotes the unchanged spatial identity multiplier of the homogeneous response, not a new constant physical history or altered initial data. No contact contribution is dropped from that anchor.
