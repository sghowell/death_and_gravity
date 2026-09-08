# Primary-source audit and withheld expansion

The accepted geometric method is supported by
[Garcia-Recio and Salcedo, arXiv:1904.06154](https://arxiv.org/html/1904.06154),
in particular their scalar determinant construction (23)--(27).
Their dimensional parameter is (d-4)/2; the executable formulas
convert to (4-d)/2. Their quadratic two-derivative expansion
(114) supplies an independent check. Derivative indices precede
the two mass-tensor indices in their shorthand.

[Ruf and Steinwachs, arXiv:1806.00485v2](https://arxiv.org/abs/1806.00485v2)
gives the generalized Proca calculation, the scalar-mass
conformal check (93), and a quadratic expansion (115).
The literal 48-term fourth-order transcription in invariants.py
was checked against the original arXiv TeX source, not only the
HTML conversion. In particular its coefficient -80 is present
in the TeX. The zero- and two-derivative parts pass independent
checks, as does the flat limit of that fourth-order transcription.

## Failed control, not an accepted formula

For the literal fourth-order implementation, substitute
Y=m^2 f g and subtract the independently derived conformal Proca
result. Its weighted Euler derivative is nonzero. At H=1,
all higher H jets zero, k=0, f'=1 and all other f jets zero,
the exact residual is 237/20. These local jets suffice to
detect a nonzero differential polynomial; smooth compactly
supported test profiles can realize them near a point.

The implementation is therefore withheld and never called by
kernel.pole(4). The accepted result is instead derived from
the complete auxiliary-metric scalar heat coefficient and
passes the conformal, constant-anisotropic and frozen flat
checks. No coefficients were fitted to make those checks pass.

This audit does not assert an erratum or identify the ultimate
cause of the discrepancy in the published expansion. The later
paper confirms agreement at the general bimetric level and
displays a quadratic comparison only through two derivatives.
That does not repair or certify this literal fourth-order
transcription. Its failing residual remains in the report as
a reproducible control, with its scope explicit.

The full geometric curvature basis is retained for any future
dimensional continuation. Neither paper supplies this project's
finite selected-state retarded response or its V/G/B obligations.
