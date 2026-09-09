# Exact classical nearby wavepackets

This follows the actual nearby [cone proof](../README.md).
It controls finite-frequency solutions of the same classical
linearized equations in a local linear relational matter observable.

The proof keeps the complete finite-frequency Hamiltonian and uses
a symplectic mode normalization whose leading diagonal transport
vanishes exactly. Uniform time-analytic bounds then control mode
mixing and a near-identity diagonalization. The resulting packet
estimate includes three-dimensional diffraction, nonzero sinc tails
and normalization of the exact reconstructed field.

Read:

- [Exact canonical and mode normalization](notes/canonical.md).
- [Complex actual-solution and derivative bounds](notes/analytic.md).
- [Local relational observable and reconstruction](notes/observable.md).
- [Normal-form error and three-dimensional packets](notes/packets.md).
- [Limits and source comparison](notes/scope.md).

With every P8 source root on PYTHONPATH:

    python -m p8_proca_nearby_packets.verify --check

This is read-only; omit --check to print a fresh report. The report
pins every local source and proof and fully rebuilds S6.89.
Native scientific arithmetic is unchanged. Written analytic,
Fourier and localization arguments have exact executable anchors;
they are not proof-assistant formalization.
