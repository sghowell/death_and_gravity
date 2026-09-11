# Independent calculations and primary-source cross-checks

The independent tests compare the logarithmic threshold
integral with its real bubble formula on both sides of
s=4. They compare boundary-subtracted physical triangle
and box values with analytic radial expressions and,
in the Euclidean region, with direct positive parameter
quadrature. The tests retain the threshold imaginary
part and detect early absolute-value and missing-factor
shortcuts. Explicit four-by-four Dirac matrices test the
resolvent-product norm and the active-flavor normalization.

For a box the independent radial substitution H=r/(1+r)
gives Delta=(delta+Mr+C r^2)/(1+r)^2,
C=M-z y(1-y), and measure r dr. With
d=sqrt(M^2-4C delta) and
L=log|(M+d)/(M-d)|+i pi 1_(delta<0),

    integral_0^infinity r dr/(delta+Mr+C r^2-i0)^2
       =(M L-2d)/d^3.

This provides a check of the physical boundary value,
not an absolute integral through its real pole. Euclidean
positive-delta checks independently verify normalization.

Primary sources were checked on 2026-09-10:

- E. T. Tomboulis, [Causality and Unitarity via the Tree-Loop
  Duality Relation, arXiv:1701.07052](https://arxiv.org/pdf/1701.07052),
  sections 4.1 and 5, derives local Hermitian-vertex cutting
  identities and discusses their multiloop extension. The
  positive-physical-state restriction in gauge theories is
  a separate issue; no cut gauge state enters our order-two
  window.
- C. de Rham, S. Melville and A. J. Tolley,
  [Improved Positivity Bounds and Massive Gravity,
  arXiv:1710.09611v2](https://arxiv.org/html/1710.09611v2),
  sections 2.1-2.3, equations (12)-(16) and (25)-(31),
  distinguishes the exact dispersive relation from its
  finite-order cut approximation.

These sources cross-check normalization and scope. No
Galileon/massive-gravity power counting or high-energy
bound is imported. The local canonical parent, not an
arbitrary nonlocal action, supplies the unitarity
calculation; S6.160 supplies only its matched physical
source representation.

Tests and written estimates are not proof-assistant
formalization or independent peer review. Native reports
pin the proofs as well as the algebra and test sources.
