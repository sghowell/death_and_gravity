# A finite first spectral outer-bubble contribution

P8 S6.135 derives the positive spectral representation of the first
on-shell fermion propagator insertion, then performs the local-quartic
outer-bubble bound with its proper inner references and one overall
constant subtraction. Both light-line positions and all three channels
are retained.

The literal family's positive forward second coefficient is less than
10^-1418 and less than 10^-819 times the tree reference coefficient.
This is not the complete two-loop amplitude or original P8 closure.

See FORMULATION.md and all six proof notes. Read-only native replay:

    .venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_fermion_spectral_bubble.verify

The report recursively rebuilds the frozen S6.134 ancestor and compares
all report fields and source hashes. Runtime allowances do not change
scientific SymPy operations; only the separately documented full regression
uses its audited exact GCD adapter.
