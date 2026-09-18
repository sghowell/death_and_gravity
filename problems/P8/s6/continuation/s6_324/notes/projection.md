# All proper probability-overlap terms

The eight anchored components gS are the Boolean Mobius transform
of the eight compatible face values. The zeta and Mobius matrices
are exact inverses. Each gS depends only on energies in S and
vanishes whenever any of these energies is set to zero.
Consequently R(gS*conj(gT)) is the same product if S union T is
the full three-label set, and is zero otherwise.

For each label the union-covering ordered pair has three choices:
membership in S only, T only, or both. Hence there are3^3=27
ordered products in R(|G3|^2). Fifteen contain g123 on at least
one side and sum to |G3|^2-|H|^2 with H=P G3. The remaining12 are
 T=2 Re[ga*conj(gbc)+gb*conj(gac)+gc*conj(gab)
        +gab*conj(gac)+gab*conj(gbc)+gac*conj(gbc)].
Equivalently T=|H|^2-P(|G3|^2), and exactly
 R(|G3|^2)=(|G3|^2-|P G3|^2)+T.

The module checks the full identity for16 independent real and
imaginary component coefficients, the inverse matrices, all three
face annihilations and the12/15 partition. Tests independently
apply inclusion-exclusion to a complex polynomial and verify all
six permutations of the labels. A sample with only gab=gac=1
gives T=2; the singleton/pair terms alone would give zero.
Thus the pair/pair terms cannot be omitted.

As a lower-order check, for two faces X,Y and corner V,
 |X+Y-V|^2-|X|^2-|Y|^2+|V|^2
 =2 Re[(X-V)*conj(Y-V)].
Amplitude-baseline and probability-baseline subtraction are
different already at two real emissions. Their finite transfer
requires quantitative face bounds; it is not a formal equality
of the two prescriptions.

Finally a generic nonnegative squared amplitude does not make
its rectangle nonnegative: R((1+ab-c)^2)=-2abc. This example
guards the logical scope only; it is not asserted to be the
original physical amplitude. The certified measure stays signed.
