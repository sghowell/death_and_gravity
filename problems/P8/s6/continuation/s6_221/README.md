# P8 S6.221: tame classical scalar propagator

The complete current coefficient-sector scalar system has a uniform polynomial-momentum comparison through the unit bounce window: exp(1e29)(1+|P|²)^6. The two-chart proof retains full finite-transfer retuning and weighted time boundaries. It loses twelve spatial derivatives; it does not establish the coupled quantum inverse or original P8 closure.

Read [FORMULATION](FORMULATION.md), [charts](notes/charts.md), [energy](notes/energy.md) and [comparison](notes/comparison.md). The seven proof notes and independent tests are source-hashed. The report is read-only.

Ordinary replay: `.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_221`.
CLI replay: `.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_scalar_tame_propagator.verify`.
Full replay retains the frozen S219 helper-path allowance:
`PYTHONPATH=problems/P8/s6/continuation/s6_219/tests .venv/bin/python -u scripts/p8_replay.py full`.
Only full regression uses the separately audited exact-GCD adapter.
