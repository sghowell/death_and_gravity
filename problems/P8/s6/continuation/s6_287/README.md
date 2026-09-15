# S6.287: ordinary gravity Ward bridge

The complete on-shell zero-transfer gauge-fixing difference vanishes.
All-topology parameter bounds then justify the ordinary pure-GR
one-loop F1 pole/finite zero-transfer limit. The proper endpoint and
scalar LSZ legs cancel the entire raw residue derivative.

At nonzero transfer the explicit remaining IR pole is still present.
This is a regulated charge-normalization result, not a physical
IR/detector theorem or full P8 closure.

Start with [FORMULATION.md](FORMULATION.md), then
[whole gauge variation](notes/gauge.md),
[all-domain continuity](notes/continuity.md) and
[raw sign and charge cancellation](notes/soft.md).
[Matching](notes/matching.md) and [scope](notes/scope.md) retain every
unresolved physical obligation.

From the repository root, independent checks use:

```sh
.venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_287
.venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_ordinary_gravity_ward_bridge.verify
```

Those modes and the native report use original SymPy GCD. Only a complete
captured FULL replay may use the separately audited exact adapter.
The raw certificate is read-only; validation history and complete-FULL
coverage are external to this frozen manifest. Original V/G/B/P8 remain
OPEN, with scoped P8(a) and all historical qualifications unchanged.
