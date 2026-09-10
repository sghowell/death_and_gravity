# P8 S6.112 — actual one-loop light pole

This child fixes the physical light pole and residue of the **new polynomial
vacuum model** already defined in S6.110. It does not transfer the earlier
cosmological state, counterterms or kinetic parent.

The full reduced Hessian yields the mixed light-heavy two-point bubble,
including the correct one-loop half-trace normalization. Local tadpole and
affine counterterm pieces cancel under one fixed on-shell subtraction.
The truncated one-loop inverse has pole mass one, residue one and no extra
zero in the complex unit disc about that pole.

The required finite local kinetic subtraction has absolute value below
3 times 10^-208. The corresponding once-only conversion of the potential
scheme changes its quadratic curvature from one to 1-p, with
0<p<10^-405. The fixed quartic subtraction and S6.110 positive one-loop
Phi^6 bound remain valid to this order.

See [the formulation](FORMULATION.md), [normalization](notes/normalization.md),
[analytic estimates](notes/analytic.md), [subtraction conversion](notes/subtraction.md)
and [scope](notes/scope.md). Run the repository's ordinary pytest or
`python -m p8_vacuum_light_pole.verify --check` with its P8 source roots.
The replay is read-only and recursively rebuilds the frozen ancestry.

Original P8 and the complete V/G/B contract remain OPEN.
