# P8 S6.216: corrected retarded phase and ordered scalar UV

This checkpoint diagnoses the inherited creation/annihilation phase mismatch and supplies the corrected full spatial UV coefficients and original local finite difference, including all three ordered scalar channels.

It does not restore S213's physical current assembly or S215's known tracefree input, and does not complete the scalar current, homogeneous trace anchor, reduced inverse or original P8.

Start with [FORMULATION.md](FORMULATION.md), [the phase proof](notes/phase.md) and [the scope](notes/scope.md). The two independent test helpers are included in the frozen source manifest. Every historical source/test/report byte is preserved.

Replay with the repository's original-SymPy ordinary/CLI entry points. The separate full regression alone uses the audited exact-GCD adapter. Original V/G/B and P8 remain OPEN.
