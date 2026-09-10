# P8 S6.113 — same-model one-loop forward coefficient

In the new S6.110 polynomial vacuum, the complete one-loop correction to
the mass-one forward scalar coefficient obeys the exact rational bound

|b2_one_loop| < 10^-6 b2_tree,   b2_tree=4 lambda>0.

Thus the tree-plus-one-loop coefficient is strictly positive in the
explicit fixed subtraction scheme. This is not an all-orders positivity
verdict, a UV completion, a finite-gravity contour bound or bounce matching.

The calculation uses all three light-loop channels with the exact nonlocal
quartic vertices obtained by integrating the Gaussian heavy field.
This retains the original heavy/light bubble, triangle and box diagrams.
Only external shifts are expanded uniformly; radial loop momentum still
runs from zero to infinity. The angular remainder is controlled by a
complex-disc Cauchy estimate, and the radial part by a subtracted
second-derivative bound with ultraviolet cancellation made explicit.

The local full-model counterterms implement one declared subtraction and
retain the parent's fixed constant-field quartic condition. The light pole
and residue are those fixed independently in S6.112; no quantum transfer
through the classical S6.111 field map is used.

See [formulation](FORMULATION.md), [normalization](notes/normalization.md),
[analytic domain](notes/analytic.md), [angular bound](notes/angular.md),
[radial bound](notes/radial.md), [subtractions](notes/subtraction.md) and
[scope](notes/scope.md). Native read-only replay is
`python -m p8_vacuum_forward_loop.verify --check` with repository source roots.

Original P8 remains OPEN.
