# S6.242 — full prepared-heavy-state selection response difference

The unchanged S240 SLE now has a complete state-selection response comparison against its exact W6-initialized KG basis. All lapse, shift and noncommuting spatial-metric contacts are retained. An independent covariance calculation checks the full ordered quantum kernel, including the reflected reverse term in odd shift cross channels.

The ENTIRE response difference, together with its algebraic summand of the existing fixed reference profile, is below10^-1280 after kappa0 normalization in the specified physical ADM and common-clock weak graphs. All internal and external momenta are included; no time derivative of the external directions is required. The full two-leg frequency-cutoff error is below10^-1180/K on the same clock graph, including the nonzero finite-cutoff chart contact.

This is not the exact-comparison state's full renormalized response or a quantum inverse. The physical state, finite prescription and full profile sum are unchanged. Interacting loops, nonlinear bounce, quantum gravitational limit, UV, Regge and original P8 remain open.

See [formulation](FORMULATION.md), [Hamiltonian](notes/hamiltonian.md), [ordered kernel](notes/ordering.md), [state comparison](notes/state.md), [full bounds](notes/bounds.md), [profile and clock contact](notes/profile.md), [scope](notes/scope.md) and [validation](notes/validation.md).

Read-only checks:

```sh
.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_242
.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_heavy_state_response.verify
PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full
```

Only full regression uses the audited exact-GCD adapter. Native/direct/ordinary/CLI retain original SymPy. The algebra and written proofs are not FORMALIZED.
