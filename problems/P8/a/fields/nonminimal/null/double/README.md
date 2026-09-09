# P8-A.22: physical double-null smearing

The [formulation](FORMULATION.md) gives a finite nonoptimal physical
null-stress difference inequality on a timelike two-plane, including
a one-sided Wick upper cap. An explicit normalized compact sampler
has a fully reconstructed coefficient. The conformal FLRW transport
keeps physical volume, affine normalization and the actual scalar
reference. No null focusing theorem or original P8 closure is claimed.

Proofs: [plane and positivity](notes/plane.md), [boost and compact
profile](notes/boost.md), [curved transport](notes/curved.md),
[scope](notes/scope.md), and [primary-source context](notes/sources.md).
The read-only report is certificates/physical-double-null-qei.json;
its verifier is src/p8a_double_null/verify.py.

Native checks leave scientific SymPy unchanged. Every report field
has omission/replacement controls. The complete P8 snapshot regression
uses the separately audited exact GCD runner. A.21 and its physical
scalar ancestry are replayed, not used as a null-line inequality.
