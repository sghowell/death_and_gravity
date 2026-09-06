# P8 continuation audit after A.11 and S6.6.COMPOSITE

Recorded 2026-09-06. This updates, without rewriting, the
[preceding audit at the A.10/S6.5 checkpoint](assessment-2026-09-06-p8-closure-audit.md).
It is an assessment, not an amendment of the original problem or a
new theorem. Original P8 remains open.

## What changed

| Track | New result | Important boundary |
|---|---|---|
| [A.11](../problems/P8/a/applicability/reference/asymptotics/backreaction/response/remainder/residual/existence/preparation/FORMULATION.md) | A complete actual-state SEE map, exact constraint-preserving preparation, strict contraction and smoothness on one specified short slab | Only the final half-slab is source-free; neither whole-target continuation nor the old QSEI's transfer is proved |
| [S6.6.COMPOSITE](../problems/P8/s6/matching/composite/FORMULATION.md) | A regular free-canonical composite bounce, and exact local CD with a reconstructed potential on `abs(T)<=tau/64` | The physical metric/shared coupling differ from S6.5; background existence does not prove mode stability or controlled matching |

The certificate hashes at this checkpoint are:

    A.11: 8dba8215c5a28e2a6379c09e4d452838449065c396e10888aa733b186799456f
    S6.6.COMPOSITE: 7708eebb7f10be8bad523029aa7a53921c943fd6f7df66691b71651d090ec06f

Their original ancestor subtrees remain unchanged. The source-hashed
proofs distinguish exact symbolic/rational checks from the written
analytic existence, smoothness and ODE arguments; none is Lean-formalized.

Validation: **1,210 P8 tests passed** in 405.35 seconds, including 75
A.11 tests and 39 S6.6 tests. Both read-only certificate replays and the
full P8 Ruff check pass. A.11 hashes 19 files and has 18 exact symbolic
identities plus its independent Fraction replay; S6.6 hashes 16 files
and has 75 symbolic and 11 independent coefficient identities. Unrelated
P4/P9 work is excluded from this checkpoint.

## P8(a): local actual applicability is now established

The former gap was not just a missing numerical coefficient. An actual
energy constraint had to hold, and finite-regularity existence could not
be promoted to a smooth Hadamard solution. A.11 addresses both by keeping
the original past and co-evolving the metric/state during a smooth,
exactly conserved preparation.

The source's weighted density is prescribed from the old actual defect.
Its pressure is fixed by conservation on the **unknown** metric. Integration
by parts removes cutoff derivatives from the uniform map bounds without
dropping that pressure. The full map retains the gravitational
`1/(60delta)` coefficient, anomaly, fixed subtraction prescription,
all old mode history and the positive inverse pole. Every differentiated
Picard map has the same strictly contractive highest-jet block, which
supplies smoothness on the same interval.

This is genuine progress beyond a small residual or a conditional Banach
test. Its cost is a very short dimensionless interval and a preparation
stress which need not obey energy conditions. It is not a globally
compact source on the unchanged off-shell radiation past. No bound for
all higher pressure/curvature jets is claimed.

The next concrete obligations are quantitative two-point/QSEI control
on an exact solution, a domain long enough for nonvacuous focusing, and
the original realistic-field requirement. A.9 cannot simply be relabelled
as a theorem on the new metric. Repeated local existence alone does not
control accumulated response or overcome the inverse's positive pole.

## P8(b): a different parent class survives the background screen

The old HR exclusion assumed separately conserved NEC sources and
physical metric `g`. The new action couples one canonical scalar to the
positive-root composite metric `g_e`. The induced sources of `g` and `f`
are not separately conserved. The old polynomial root theorem therefore
does not apply to its pressure-dependent branch factor.

Both lapse and scale equations, both Bianchi balances and branch
intersections are retained. The free example crosses a double pressure
root smoothly. An independent fourth-jet invariant shows it is not CD
even after changing the proposed CD time scale. The separate CD
construction reconstructs a local scalar potential; it is not a solution
for frozen free M1 matter or a controlled match to the old CD/M1 action.

The next screen is the actual rolling quadratic theory, retaining tensor,
vector and both scalar channels and their constraints at the bounce.
A vacuum mass, a minisuperspace constraint or a published high-frequency
kinetic sign does not establish the full rolling gradient and cutoff
conditions. Composite coupling also has a known generic constraint/EFT
boundary; the [source audit](../problems/P8/s6/matching/composite/notes/sources.md)
keeps it separate from the background construction and from corrected
trimetric actions.

An independent parent interaction scale must not silently be identified
with the CD duration when testing a broader family. Conversely, modifying
a frozen candidate to pass one sign check creates a new candidate whose
remaining equations and matching conditions still require verification.

## Closure decision

The frozen 32-row linear-principal classification and these two new
bounded gates are complete as stated. The realistic-field focusing and
S6 controlled-matching/positivity obligations remain. Neither a full
all-orders UV completion nor arbitrary nonlinear stability should be
added as an unannounced replacement for the original contract.

No missing credential, external approval or user preference has emerged.
The next work is mathematical, with falsifiable checks on both tracks.
Original P8 is **not finished or closed**.
