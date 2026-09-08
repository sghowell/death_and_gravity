# S6.50: physical vector energy and an all-momentum evolution bound

See [FORMULATION.md](FORMULATION.md), the
[physical energy derivation](notes/energy.md) and the
[exact mode comparison proof](notes/comparison.md).

The certificate concerns the unchanged retained vector on the actual
clock, with explicit fourth-order Gaussian initial data on a compact
time interval. It does not establish full renormalized energy,
all-order Hadamard admissibility, V/G/B or original P8 closure.

The read-only entry point is `p8_vector_state.verify --check` with
the P8 source roots on PYTHONPATH. Replay fully rebuilds the frozen
parent and checks every local source hash. Ordinary tests do not use
the broad exact GCD adapter. Repository-wide regression uses the
established independently controlled runner without frozen exclusions.
