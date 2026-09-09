# P8-A.23: a conditional physical-QEI null theorem

The [formulation](FORMULATION.md) gives global FLRW null affine
incompleteness within 2tau using the actual conformal scalar
double-null QEI, explicit finite-plane state and geometric caps,
and the actual semiclassical equation. No homogeneous-state or
pointwise NEC/SEC assumption is inserted. Original P8 stays open.

Proofs cover the [physical source reduction](notes/plane.md),
[outgoing index identity](notes/index.md), [full constants](notes/constants.md),
[complete geometry control](notes/control.md), [actual thermal
past](notes/thermal.md), and [scope](notes/scope.md).

The read-only certificate is
certificates/double-null-qei-flrw-incompleteness.json and its verifier
is src/p8a_null_focusing/verify.py. Native scientific runs preserve
SymPy. The complete snapshot regression separately self-checks its
exact GCD adapter. The report rebuilds A.22 and its full ancestry,
pins all source/proof/test files, and rejects mutation of every field.
