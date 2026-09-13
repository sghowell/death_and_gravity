# S6.245: the complete actual heavy-scalar spatial response

This continuation establishes the entire first spatial-metric response of the unchanged QG2 heavy SLE on its reference FLRW slab, at every external Fourier transfer. Both trace and tracefree directions, the complete state error, the sixth time remainder, every UV-subtracted endpoint and the unchanged finite/profile terms remain.

For prepared sources, define
\[
 Z_{136}(G)^2=\sum_{j=0}^{13}\|(1-\Delta)^3\partial_t^jG\|_{L^2(I\times\mathbb R^3;F)}^2,
 \qquad I=[-1/2,1/2].
\]
The complete normalized response satisfies
\[
 |\mathcal R(D,G)|/\kappa<10^{-350}\|D\|_{L^2}Z_{136}(G).
\]
Its restored, UV-subtracted original-frequency two-leg projection satisfies
\[
 |(\mathcal R-\mathcal R_K)(D,G)|/\kappa
 <10^{-250}K^{-1}\|D\|_{L^2}Z_{136}(G),\qquad K\ge2\sqrt n.
\]
The detector needs no spatial derivatives. No finite external-transfer ball is imposed.

The projection is defined by the complete subtracted decomposition in [assembly](notes/assembly.md), with all unaveraged UV terms retained before masking and all finite local/profile terms restored. It is not the bare cutoff Gaussian current or a finite-cutoff Ward identity.

Read [FORMULATION](FORMULATION.md) for the exact claim, [domain](notes/domain.md), [state](notes/state.md), [time](notes/time.md), [endpoints](notes/endpoints.md) and [dimension](notes/dimension.md) for the estimates, and [validation](notes/validation.md) for the independent checks.

The native certificate is read-only. Its 17-source manifest and 20 report fields are replayed without changing the frozen inputs. Native/direct/ordinary/CLI use original SymPy; only full regression uses the separately audited exact-GCD adapter.

This is not the full ADM/common-clock assembly, a quantum inverse, nonlinear inhomogeneous feedback, interacting-loop control, a physical UV/Regge completion or original P8 closure.
