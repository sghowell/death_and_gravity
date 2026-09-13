# S6.253: homogeneous physical-background vertices

The same full parent now has an explicit jointly reduced scalar and
longitudinal-Proca Gaussian for homogeneous background probes. Keeping the
source square exposes mixed first vertices and second contacts absent from
a sum of independent source-free light and vector determinants.

The complete canonical6x6 first/second matrices have an exact18-row physical
substitution map. All remaining physical modes and their existing state
maps are retained. Aligning the temporal vector with its classical source
is distinguished from holding W0 fixed: their second variations differ by
the full normal-vector onepoint contact. The original P8 frontier is not
closed by these finite Gaussian results.

See [FORMULATION.md](FORMULATION.md) for the precise scope and written proofs.
From the repository root, replay with:

```sh
.venv/bin/python -u scripts/p8_replay.py ordinary problems/P8/s6/continuation/s6_253
.venv/bin/python -u scripts/p8_replay.py cli p8_vacuum_affine_physical_background_vertices.verify
```

The certificate is read-only. Scientific sources, tests, proof notes and
native report are frozen together; later validation results belong in the
external dated assessment. No older scientific file is modified.
