# Pinned sources and physical boundary

* A.11 `smooth-preparation.json`, SHA256
  `8dba8215c5a28e2a6379c09e4d452838449065c396e10888aa733b186799456f`,
  supplies the exact unchanged map, source, state, initial constraint,
  original uniqueness ball, A.8 plateau bounds and smooth flat-start
  Dyson/Hadamard proof. Its source duration is read as L0, not replaced
  by the new existence interval L.
* A.13 `causal-resolvent.json`, SHA256
  `1569abacc14b4c53c951fc10e07eea7e5f17625812705bfb5709772a26261d97`,
  supplies the exact full-map preconditioner and the complete causal
  inverse, including the growing pole, with weighted norm <=240/769
  at sigma=2*10^6. A.14 does not delete that pole or treat this inverse
  alone as the extension theorem.

The new theorem consists of the complete weighted remainder estimates,
pointwise barriers, self-map, smoothness adaptation and same-source
uniqueness comparison proved here. The actual infinite-frequency state
estimates are the pinned A.10/A.11 lemmas, not a finite mode calculation.
No new finite counterterm, order reduction or state choice is made.

For context, [Meda--Pinamonti--Siemssen, arXiv:2007.14665](https://arxiv.org/pdf/2007.14665),
Section 5, motivates inversion of the retarded logarithmic term before
using a fixed point. Theorem 5.9 and Remark 5.2 do not supply this new
quantitative massless smooth interval. The present proof uses the
already audited actual-state lemmas and derives its own full constants.
As emphasized in the A.13 primary-source audit, a weighted estimate does
not remove runaway poles or authorize a different causal prescription.

The cutoff's physical high-derivative validity and quantum stress
fluctuations are not assessed. A longer mathematical solution is not
automatically a cosmological-scale application. The next required
calculation is a fresh all-Hadamard QSEI/reference bound and a focusing
duration check on this actual extended geometry.
