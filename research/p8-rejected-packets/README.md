# Rejected P8 packets: immutable evidence, not accepted certificates

The s6_279 directory is the byte-for-byte archive of an UNPUBLISHED
packet rejected by independent native validation on2026-09-15.
It was never part of an accepted complete P8 regression or public claim.
No source, test or native-report byte has been repaired in place.

Its private and fresh source preflights each produced130 named checks
and391 scalar entries; all321 science tests passed. The independent
full-ancestry native run produced139 named checks and409 scalar entries.
The native raw report is132764 ASCII bytes, SHA256
a45d6d297e8ecea8f6049983b69c8914125ef2fde89cea537fe578d34064d256.

Cause: S278 source.full_vacuum returns a cached mutable check dictionary.
S278 source.data extends that same dictionary in place by9 checks and18
scalar entries. S279 copied every entry, so its own report depended on
whether the parent's full data packet had been built first.

All18 source hashes matched the accepted private preflight, and native
raw-byte/ASCII/field/source checks passed before the count comparison
failed. This is a report reproducibility defect, not a nonzero algebraic
residual. The archived tests deliberately retain their original contract
and are not advertised as passing native certification.

The rejected packet has been moved outside the accepted P8 replay tree,
with all18 scientific sources and the raw native report preserved.
This does not remove any previously accepted P8 certificate or test.
The completed809-file S277/S278 snapshot excluded this later packet,
as does the completed813-file,71,032-test S280/S281 snapshot.

The corrected [S280 successor](../../problems/P8/s6/continuation/s6_280/README.md)
selects the intended original ten-check source set explicitly and copies
it without aliasing. Cold, forced-parent-warm and fresh-repository
preflights agree exactly. A new untouched original-SymPy native process
first reproduced S278's84449-byte report exactly, then returned S280's
expected130 named checks and391 scalar entries. Independent ordinary,
CLI and complete-regression gates passed. The rejected evidence above
remains unchanged; the successor does not certify the archived packet.

The old native process was retired with its original protected snapshot.
No helper, frozen source or raw report was repaired in place.
