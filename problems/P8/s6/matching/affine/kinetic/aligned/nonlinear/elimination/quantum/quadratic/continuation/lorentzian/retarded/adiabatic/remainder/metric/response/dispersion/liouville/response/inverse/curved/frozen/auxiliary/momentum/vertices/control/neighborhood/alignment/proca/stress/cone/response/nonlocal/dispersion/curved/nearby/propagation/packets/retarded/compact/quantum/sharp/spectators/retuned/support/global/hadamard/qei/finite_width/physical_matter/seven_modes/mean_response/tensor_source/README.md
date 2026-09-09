# Actual tensor clock-source and global mean-response certificate

The new positive TT state difference retains the actual off-clock
tensor coupling, induced scalar matter, the nonlinear center
spatial chart and the clock-source Ward identity. Global bounds
give a complete bouncing first-order metric representative,
explicitly not an exact quantum solution.

See [formulation](FORMULATION.md), [action](notes/action.md),
[state](notes/state.md), [chart](notes/chart.md), [mean](notes/mean.md),
[clock exchange](notes/clock.md), [bounds](notes/bounds.md) and
[scope](notes/scope.md).

Scientific replay uses the repository Python environment, all
P8 src directories on sys.path and

    python -m p8_tensor_mean_source.verify --check

The certificate CLI is read-only. Native science and ordinary
replay do not monkeypatch SymPy. The separate complete P8 runner

    .venv/bin/python -u scripts/p8_snapshot_regression.py

captures the test-file snapshot and uses the separately audited
exact GCD regression adapter. That adapter is not a substitute
for the native scientific/report replay. Source hashes include
this formulation, every note, own Python source and test file.

The saved report is
[actual-tensor-clock-mean-response.json](certificates/actual-tensor-clock-mean-response.json).
Original P8 remains open.
