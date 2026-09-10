# P8: controlled first gauge light-cut subtraction

Original P8 remains OPEN. This follows the
[finite-mass first-cut audit](assessment-2026-09-10-p8-finite-mass-first-gauge-cut.md).
No user intervention is needed for the current research.

## A specific local repair

[S6.130](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/quadratic/continuation/lorentzian/retarded/adiabatic/remainder/metric/response/dispersion/liouville/response/inverse/curved/frozen/auxiliary/momentum/vertices/control/neighborhood/alignment/proca/stress/cone/response/nonlocal/dispersion/curved/nearby/propagation/packets/retarded/compact/quantum/sharp/spectators/retuned/support/global/hadamard/qei/finite_width/physical_matter/seven_modes/mean_response/tensor_source/scalar_source/global_response/vacuum/affine_domain/renormalizable/offshell/light_pole/forward_loop/elastic_cut/spin2/spectral/running/quantum_map/euclidean/two_loop/finite/insertions/double_bubble/wineglass/finite_contact/light_pole/canonical/uv_screen/finite_mass/subtraction/FORMULATION.md)
keeps the prospective GY14-unbroken candidate and its
first two-gauge channel at three parent loops. It
subtracts both known low cuts explicitly, instead of
discarding a tiny but nonzero discontinuity or importing
the old scalar gap hypothesis.

The transition remainder is bounded on the complex
disc |s-2|<=5. Scalar Bose symmetry and Lorentz
invariance let the continued invariant amplitude pass
through the square-root coordinate thresholds. The
Ward identities give the soft endpoint factors. For
z=a y^2/(Q mF^2), epsilon=4 times 10^6/mF, and
K=dA z^2/(9pi^2), the holomorphic spectral error satisfies

    |delta rho| <= 9pi K eta,
    eta=(5/2)epsilon+epsilon^2/8.

This is complex-domain control, not a real-axis estimate
used to justify an unbounded derivative.

Define the known low-cut removal, initially off the real
axis, by

    F_low(s)=(1/pi) integral_0^6 rho(x)
               [1/(x-s)+1/(x-(4-s))] dx.

Its crossed boundary signs remove exactly the specified
first-channel discontinuity locally. A removable
divided-difference representation controls each analytic
boundary germ. The two germs can differ by a crossing-odd
term; their even center derivatives agree. This does not
restore an all-loop mass gap.

## The coefficient and what it does not determine

The leading center second coefficient is

    b2_low = K[2 log(2)-21/4] + delta b2_low.

The explicit complex bound, two finite integrals and
logarithmic endpoint terms give |delta b2_low|<=288K eta
by Cauchy's estimate on the unit disc. No singular
principal-value integral is differentiated naively.

At the named mass mF=10^200, the relative coefficient
error is below 10^-190, and rational logarithm enclosures
prove

    -4K < b2_low < -3K,
    |b2_low| < 4 times 10^-1620.

These are bounds on the subtraction coefficient, not
the full new-model forward amplitude. Adding an analytic
crossing-even quadratic term would leave the specified
cuts unchanged while changing the center coefficient.
That diagnostic demonstrates missing matching data; it
does not authorize retuning the actual parent action,
and it would affect a separate high-energy growth test.

The regular matched amplitude, later light cuts and
high-energy contour therefore remain distinct inputs.
The calculation excludes neither the candidate nor a
whole ladder row.

## Independent verification

The native report pins 18 source/proof/test files and
20 fields: 80 named exact identities and scalar
entries, 29 proof gates, 9 controls and 24 rejected inputs.

Final private science: 164 tests in 20.82 seconds.
Independent repository science: 164 tests in 20.93 seconds.
Ordinary read-only replay: 189 tests in 2078.71 seconds.
Independent native CLI: passed.
Complete P8 snapshot: 18173 tests in 2867.37 seconds.

All 513 captured test files are present and unchanged.
Path-list SHA-256:

    6e8dced4ee4bc0bc4b603e671056d2c955879369139f37bad010440cbb1f2818

The exact GCD adapter passes 128 original tuple comparisons.
Full-run counters: 33975 domain fallbacks, 7340 exact
descents, 4 mixed fallbacks. The snapshot predates
S6.131's two test files. Native, ordinary and CLI
scientific processes use unmodified SymPy; the adapter
is confined to independently audited full regression.
Native ancestry reuse follows immutable source/report
hash comparisons.

The 34763-character native report was transferred
losslessly and all 18 source hashes were verified.
Report SHA-256:

    0698704470509fa022fa4ccf098a27c20b4724b1235eb7387785fc735ce7bb78

Exact 22-file staging includes only the checkpoint,
this root audit, CLAIMS.md and README.md. It verifies
all staged native source/report bytes and excludes
nested children and unrelated P4/P9 work.
Complex-domain, boundary-gluing and Cauchy arguments are
written proofs, not proof-assistant formalized or peer reviewed.

## Active continuation and original obligations

S6.131 is frozen and native-certified, with 202 tests
passing independently in private and repository science.
Its fresh ordinary, CLI and complete-snapshot replays
remain in progress. It derives the complete one-loop
fermion local increment, retains the large mass and
vacuum-energy references, and makes the canonical field
map explicit. Its positive finite-domain potential is
only the scalar tree plus one-loop fermion functional,
not the full new-model quantum potential.

The next matching task is the momentum-dependent
four-scalar fermion contribution and its proper
reference/field dictionary. The old complete scalar
two-loop budget cannot simply be assigned to the new
candidate. Required truncation errors, justified V
contours, finite-gravity IR/Regge allowance and the
common bounce-parent field/state/cutoff domain remain.

These are the adopted finite-EFT and necessary-positivity
obligations, not a newly imposed demand for a constructed
all-orders UV theory. Scoped P8(a) and A.20-A.23 stand
unchanged. Original P8(b) and original P8 remain open.
