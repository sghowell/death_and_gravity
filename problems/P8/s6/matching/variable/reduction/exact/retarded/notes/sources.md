# Evidence and independent-check boundary

The immediate source is the frozen S6.33 exact own-f background report,
recursively rebuilt through S6.32/S6.30 to the unchanged S6.20 action.
No frozen coefficient, physical metric, canonical clock or free-chi
coupling is changed. The test geometry remains an off-shell prescribed
physical metric, not the original full rolling-g solution.

bounds.py explicitly inherits the exact-Fraction background enclosure
from S6.33. It does not claim a second independent reconstruction of that
background. It independently checks the smaller physical-domain inclusion,
the complete denominator polynomial identity, positive coefficients and
monotone rational inequalities for k and s. Its tests recompute the
polynomial convolution and integer cross-multiplied bounds.

green.py and notes/green.md derive the first-zero positive-flux argument,
two Volterra lower bounds and the smooth-pulse response estimates. The
proof is continuous in source time, observation time and the admitted
finite action parameter. No fitted trajectory or numerical Green solve
is used to assert the theorem. The explicit pulse is proved smooth at all
joins by a general derivative-polynomial recurrence, not only finite jets.

The separately authored action.py derives physical/inner time and action-
measure factors, the momentum system and the original physical-g equation
readout from the literal tensor action. It recalculates the rational Green
constants without importing green.py. Independent tests include variable-
kinetic and constant-coefficient exact kernels, the flux jump, an omitted-
drift control, exact nonunit physical scales, the factor-two variational
normalization and explicit rejection of a longer uncertified interval.
The wrapper compares coefficient/Green domain constants and independent
response constants before building the report.

This is a CERTIFIED written-proof/exact-arithmetic result, not a kernel-
checked formalization. The boundary/state prescription is part of the
theorem. Changing it can change the response; shrinking-duration metric
inputs are not fixed physical low-frequency matter sources. The causal
equation is not silently identified with a single-copy effective action,
and its original-g residual is not silently identified with a permitted
conserved matter stress. These distinctions prevent promotion to a full
CD/UV-matching or original P8 theorem.
