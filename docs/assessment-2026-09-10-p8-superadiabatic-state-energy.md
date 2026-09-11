# P8 continuation: exact-frame flat Dirac state-energy difference

Checkpoint S6.167 sharpens the complete free flat transition estimate
without changing the operator, in/out states or prior scientific bytes.
Original P8(b), V/G/B and original P8 remain OPEN.

## Result and scope

Four exact alternating rotations retain every connection term. The
analytic-disk induction gives

    |beta_p| <= 2^42 p Delta/(tau^4 (p^2+m0^2)^3),
    m0=.99m.

This is a complete quadratic-transition bound, not a truncation of its
Dyson series or a claim about all Feynman loops. At the actual parameters
the uniform amplitude is below1e-390 and the out-particle energy below1e20.

The local energy comparison retains coherence linear in beta:

    |rho_in(t)-rho_out_state(t)|
        <=8N 2^42 Delta/(3pi^2 tau^4 m0)<1e411,

for every real time, or below1e-389 relative to the named kappa scale.
The out state is not assigned zero local energy during the transition.
Identical local subtraction and finite counterterms cancel in this
state difference. Neither absolute local stress nor a bounce-density
relative error follows.

The same free flat in/out states were identified as Hadamard in S6.166.
The external exact-rotation framework is
[Barbero et al., section IV](https://arxiv.org/pdf/1805.05107);
the explicit disks, constants, full transition and coherence integrals
are derived in this checkpoint. The physical energy is the original
Dirac Hamiltonian, not a rotating-frame generator.

## Verification and recovery

Final corrected private science passed280 tests in106.43 seconds;
fresh repository science passed280 in105.39 seconds.
Ordinary replay passed305 in2199.79 seconds.
Independent native CLI replay passed.
The complete captured P8 snapshot passed27925 tests in3346.59 seconds,
final exit code0.

All587 captured test files remained present and unchanged.
Path-list SHA-256:

    51453fb35458db4c3c838497b24baa83aec5188b96ff435d473b260be5183614

The full-run exact GCD adapter passed128 original tuple comparisons;
final counters were34093 domain fallbacks,7340 exact descents and88
mixed fallbacks. Only full regression uses that adapter.
Native, direct science, ordinary and CLI use original SymPy.

The corrected native report contains20 fields and pins14 own
source/proof/test files:40 named identities,40 scalar entries,
30 proof gates,9 controls and169 rejected inputs.
Its24710 characters were transferred losslessly in three chunks.
Report SHA-256:

    01ae1949bc7d1b8eb01002e6388c2734d1bfcfb3445cb6bfcb7d7bb96512ddce

The first private/repository attempt passed its279 science tests but
failed native report encoding: a diagnostic decimal was a SymPy Float,
which the unchanged exact serializer correctly rejects. It was not an
exact-bound or proof failure, and no report, replay set or commit had
been published for that attempt.

Both frozen attempt directories were moved without changing their
bytes to the recoverable archive

    /private/tmp/p8-s767-native-encoding-failure.kYt2IR/

All28 copies were checked against the14-file manifest. A fresh private
copy changed only decimal display values to strings and added the
serialization regression check before final verification and freeze.
The exact bound was unchanged. The cold native chain rebuilt S6.166
byte-identically and then produced this corrected report. No published
scientific source or report was altered.

Exact18-file staging includes only this checkpoint, report, audit and
root ledger/link updates; unrelated P4/P9 edits and later continuations
are excluded.

## Next research

The specified absolute free flat one-loop energy and complete pressure
calculations have separate native/science evidence and fresh replay sets.
Their publication is not implied by this checkpoint.
Actual curved-state work is being developed privately.

Absolute curved and full interacting stress, physical cutoff and loop
errors, quantum target matching and controlled background response remain
open, as do vacuum contour/truncation and finite-gravity IR/Regge bounds.
Scoped P8(a) and A.20-A.23 remain closed unchanged.
