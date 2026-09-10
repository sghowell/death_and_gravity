# P8 S6.114 — actual low-energy elastic cut

The same polynomial vacuum has a strictly positive one-loop elastic
absorptive contribution on the physical interval 4<s<=6 in mass-one
invariant units. The exact nonlocal tree amplitude, not its contact
truncation, is used in the cut.

The computed low-energy contribution to the forward coefficient satisfies

lambda^2/20 < I_[4,6] < 3 lambda^2.

Combining this with S6.113's complete one-loop coefficient error gives
a strictly positive tree-plus-one-loop coefficient after subtracting
this actual low-energy cut. The conservative relative error, including
the cut subtraction, remains below one millionth of the tree coefficient.

This is a finite-order, finite-window result. It is not a full improved
positivity verdict without higher-loop and global dispersion control,
nor finite-gravity or rolling-bounce matching. Original P8 remains OPEN.

See [formulation](FORMULATION.md), [amplitude and angular integral](notes/amplitude.md),
[normalization](notes/normalization.md), [cut bounds](notes/cut.md),
[literature](notes/literature.md) and [scope](notes/scope.md).
Read-only replay is `python -m p8_vacuum_elastic_cut.verify --check`
with the repository's P8 source roots.
