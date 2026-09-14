# S6.271: explicit finite Weyl operator comparison

This continuation upgrades the S6.270 heat-symbol remainder to an actual
operator-norm estimate for the same finite, full-source regulator. It does not
identify that regulator with the original unregularized quantum theory.

The unchanged system has 48 canonical pairs, 96 phase coordinates, core radius
R = 10^20, real time |u| <= 10^-2000, the original pure covariance and exactly
the two explicitly differentiated S6.270 cutoffs. A direct normalized Gaussian
frame argument supplies the explicit Weyl operator constant 10^112; complete
mixed phase derivatives through order 196 control its heat remainder.

For each cutoff, the Weyl/calibrated-coherent volume operator difference is
less than 10^-200. The actual Weyl volume is greater than I/2 and within
10^-199 of I. The two entire propagators differ in operator norm by less than
10^-944. The two Weyl cutoff states differ by less than 10^-943 and their
volume means by less than 10^-199. The positive coherent-POVM leakage bound
for the evolved Weyl state is less than 10^-942, not zero.

These are finite local ordering and regulator comparisons. Original physical
Wilsonian matching, the unregularized interacting volume, regulator removal,
omitted loops, UV completion and nonlinear global completeness remain open.
The completed scoped P8(a) result and the S6.261 refutation are unchanged.

Read FORMULATION.md and all seven notes for the written analytic and operator
proofs. The exact certificate checks identities, rational bounds, source
bindings, ancestry and rejection controls; it does not formally verify those
written proofs. Independent numerical fixtures are diagnostics, not a
simulation of the full P8 quantum evolution.

Run the read-only original-SymPy ordinary and CLI replays with:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_271
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_weyl_operator_comparison.verify
