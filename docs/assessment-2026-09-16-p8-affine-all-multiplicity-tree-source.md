# P8 S6.315: the canonical tree source at every finite multiplicity

Original V/G/B/P8 remain OPEN. Scoped P8(a) is unchanged.

## Result and scope

This successor removes the low-multiplicity vertex truncation from the
selected four-dimensional scalar/heavy-scalar/Einstein tree source.
It derives every metric coefficient and every Einstein vertex directly
from the unchanged covariant action. For any finite labeled external
list away from internal poles, the rooted recursion enumerates each
connected tree exactly once.

The original parameters remain n=10^200/512+2, g=1/8192,
kappa=10^800 and the tuned contact C. No unknown hard loop,
evanescent coefficient or matching operator has been supplied.
The physical point API retains the original recoil, mass shells,
momentum conservation and transverse-traceless polarization rules.

## General derivation, not a finite-order extrapolation

For g_metric=eta+2*sum_i x_i*A_i/sqrt(kappa), each vertex is the
full mixed coefficient of the product of its distinct labels.
The trace-log density, inverse-metric multiplication identity,
anchored set-partition density recursion and separately labeled
Christoffel coefficients yield all matter and Einstein vertices.
The Einstein expression is the full ordered partition of the
density inverse and two connections. Linear de Donder gauge
introduces no higher gauge-fixing vertex.

Cutting the unique root edge gives the complete finite recursion.
Its exact four-scalar/N-graviton graph-count generating function is

T_N=N!*[y^N]{c*z/(2-z)^4+3*(1+b^2)*z^2/(2-z)^5},
h=y+exp(h)-1-h, z=exp(h).

Here c and b mark contact and heavy graph topologies; they are not
changed physical couplings. At c=b=1 it gives 7,47,434,5116,73444
for N=0,1,2,3,4. The analytic majorant proves
T_N<88*8^N*N! for every N>=0. Literal labeled root-cut counts
independently check the first five coefficients.

General Frobenius/Euclidean vertex bounds also follow from the
coefficient formulas. In particular, at order r>=3 the Einstein
vertex is at most r!*32^r*L^2*product_i||A_i||, times the physical
factor kappa^(1-r/2). Separate density, inverse, connection and
matter bounds retain their factorials and polarization norms.
These are written all-order proofs, not conclusions inferred
from finite coefficient tests.

## Independent calibrations and omission controls

The new formulas recover every previously implemented vertex.
A separate literal determinant/adjugate/indexed-Christoffel
calculation produces the EH5 coefficient 16438 and fifth density
coefficient -2088. Another literal adjugate calculation checks
the scalar fourth-metric vertex and all 16 inverse-density entries.

The complete new three-real coefficient contains 5116 trees.
Two exact gauge variations vanish. Removing the fifth Einstein
vertex leaves 5113 trees and a nonzero Ward defect; removing the
fourth scalar-metric vertex leaves 5110 and another nonzero defect.
The original-parameter coefficient is evaluated separately from
the diagnostic-coupling examples. No numerical zero tolerance or
modified symbolic backend enters these tests.

## What this does not establish

A graph-count bound and vertex bounds do not yet bound the
complete amplitudes near soft or collinear propagator poles.
The point API excludes internal poles, including exactly
collinear null pairs; it does not assert their continuous limits.
Later separately validated successors address additional bounds.

Neither this source nor any fixed finite tree inventory supplies
an all-multiplicity infrared-subtracted detector probability,
complete finite hard real-virtual matching, an interacting
quantum state, absolute complex Regge control, or common-parent
bounce matching. None of the original V/G/B/P8 gates is closed.

## Immutable source and native evidence

Cold 55274 passed in 247.114s and passed 657 science tests in 0.45s.
Fresh 26852 forced the complete S314 parent packets, passed in
289.238s and passed 657 science tests in 0.41s. The 18 private,
fresh and repository sources matched byte for byte before freezing.

The 18 ASCII sources total 79,211 bytes: 381 named exact checks,
531 scalar entries, 94 proof gates, 8 controls and 73 rejected inputs.
The 9 primitive records, 171 matching records and 6 historical
qualifications preserve the inherited frontier and append only the
explicitly scoped new matching result.

One native specification protected 6,354 inputs. Its 190,537-byte
ASCII report arrived in 16 chunks with 15 complete-chunk ACKs.
SHA256:
c0786edc77e68ff1c1a34b788674a782282ab96fb74be8268a345d5109840c40.
Both raw/source/count/frontier validators and the 20-field REPORT
AST passed. Native raw freeze: 2026-09-16 19:28:11 UTC.

## External acceptance

Own ordinary 5193 passed all 682 tests in 3367.15s (0:56:07).
Independent CLI 61287 returned its explicit S6.315 success message.
Both used original SymPy and exited zero. Their completion was
captured at 2026-09-17 05:54:48 UTC. Frozen source/native bytes
remained unchanged.

The complete FULL 44112 passed 88,848 tests in 6791.09s (1:53:11),
with all 889 captured files unchanged. Test-file-list snapshot SHA256:
6c99f6541e12a121814cad3c1fdf14a4f513b91c4c615c0efb0bd244364b7031.
It captured 711 namespace ancestors, 5 helpers and 128 contracts.
Counters: 68,749/7,388/94. Completion was captured at
2026-09-17 04:54:26 UTC. This snapshot includes S315 through S319.

Native/direct/ordinary/CLI use original SymPy; the exact-GCD
adapter is FULL-only. Publication requires the exact 22-file scope,
unchanged source/raw hashes and matching staged bytes, preserving
unrelated P4/P9 work. Original V/G/B/P8 remain OPEN.
