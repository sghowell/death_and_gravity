# S6.240 — specified heavy Gaussian state and fixed reference mean

The separate **CD-REG-AFFINE-ISO-QG2-H8A420** candidate specifies the added massive scalar's exact state of low energy on the full CD reference. A complete all-momentum estimate and a separately derived covariant scalar prescription bound its reference energy, pressure and five time derivatives by10^-400 against kappa0 on[-1,1].

A fixed physical scalar coefficient cancels this actual Gaussian reference stress and keeps the old Proca/M1 mean data with the S238 reconstructed affine mean. The heavy and light flat Gaussian vacuum constants are both retained. Nonconstant new vacuum jets start at degree2048, leaving S239's complete formal first-loop four-light result unchanged. The full combined clock-coefficient four-jet budget is a NEW10^-399.

These are conditional Gaussian state and mean results. The clock Hessian changes; old stability and inverse results do not transfer automatically. Full interacting light/mixed clock loops, full new response and nonlinear same-state bounce, quantum gravitational decoupling, physical UV, finite-gravity Regge and original P8 remain open.

See [the formulation](FORMULATION.md), [state](notes/state.md), [dimensional renormalization](notes/renormalization.md), [complete bounds](notes/bounds.md), [fixed profile](notes/profile.md), [scope](notes/scope.md), [primary literature](notes/literature.md) and [validation](notes/validation.md).

Read-only checks from the repository root:

```sh
.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_240
.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_heavy_curved_state.verify
PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full
```

Only the complete regression uses the audited exact-GCD adapter. Native, direct, ordinary and CLI checks retain original SymPy. Exact algebra and the written analytic proofs are not FORMALIZED.
