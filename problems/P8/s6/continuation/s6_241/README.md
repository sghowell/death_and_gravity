# S6.241 — current QG2 quadratic input and complete finite heavy local response

The existing S240 QG2 candidate now has an explicit new clock-quadratic audit. Its actual fixed combined stress satisfies every hypothesis of the generic two-chart classical comparison. The full current coefficient-sector propagator retains exp(10^29)(1+P²)^6, with twelve spatial derivatives lost. The independent added classical heavy mode has an exact energy-scaled phase bound below1000.

The complete finite heavy heat action and the entire fixed heavy-profile Hessian are derived in physical ADM variables and transferred to the common clock variables. The constant vacuum action cancels covariantly. Every lapse, shift, curvature, volume and second-metric term is retained. Both physical and same-clock weak-graph bounds are below10^-580. The complete finite tensor Euler operator has a10^-597 bound on its stated order-four graph, including the physical weighted adjoint.

These are complete LOCAL response and classical coefficient-comparison results. The remaining exact state/subtraction-dependent heavy determinant response, full quantum constraints/inverse, interacting light/mixed clock loops, nonlinear same-state bounce, quantum gravitational limit, UV, Regge and original P8 remain open.

See [the formulation](FORMULATION.md), [clock expansion](notes/clock.md), [generic comparison application](notes/propagator.md), [full local action](notes/local.md), [tensor operator](notes/tensor.md), [bounds and common-clock pullback](notes/bounds.md), [scope](notes/scope.md) and [validation](notes/validation.md).

Read-only checks:

```sh
.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_241
.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_heavy_clock_quadratic.verify
PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full
```

Only full regression uses the audited exact-GCD adapter. Native/direct/ordinary/CLI retain original SymPy. No frozen input is edited; the algebra and written proofs are not FORMALIZED.
