# Fixed ordinary Proca stress and clock profiles

S6.82; see the [formulation](FORMULATION.md),
[readout proof](notes/readouts.md), [subtraction](notes/subtraction.md),
[continuous bounds](notes/bounds.md), [new fixed profiles](notes/profiles.md)
and [scope](notes/scope.md).
The new energy readout, subtraction and finite local
matching are recomputed, without reselecting the state.

The read-only report is [fixed-ordinary-Proca-stress.json](certificates/fixed-ordinary-Proca-stress.json).
Native replay uses `python -m p8_proca_stress.verify --check`
with every frozen P8 source root on the Python path.
Original P8, full quantum response and V/G/B matching remain OPEN.
