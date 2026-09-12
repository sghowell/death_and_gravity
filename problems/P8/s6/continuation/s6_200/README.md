# P8 S6.200: full spatial flat subtraction conversion

The complete flat Proca tensor cut gives an exact conversion between
covariant and retarded time-derivative nonlocal subtraction
representatives. Three finite equal-time coefficients are spatially
nonlocal. Their uniform coefficient/tail bounds, together with the
sixth-time-derivative spectral bulk, control the entire specified flat
nonlocal representative on compact test tensors.

The physical finite local polynomial and actual curved CD contact/
endpoint matching are not determined here. Original P8 remains OPEN.

See [FORMULATION.md](FORMULATION.md), the seven proof notes and the
[read-only certificate](certificates/polynomial-vacuum-affine-flat-spatial-conversion.json).

From the repository root:

    .venv/bin/python scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_200
    .venv/bin/python scripts/p8_replay.py cli p8_vacuum_affine_flat_spatial_conversion.verify
    .venv/bin/python scripts/p8_replay.py full

Native, direct science, ordinary and CLI use original SymPy.
Only the full regression uses the separately audited exact GCD adapter.
