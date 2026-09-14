# P8 S6.271: explicit finite Weyl operator and unitary comparison

Original P8 remains **OPEN**. This is an ordering comparison for the same
finite local regulator, not original physical matching or global closure.

## Result

For the unchanged S6.270 full nonlinear source, 48 canonical pairs,
96 phase coordinates, original pure seed and two explicit cutoffs with
core R = 10^20, the complete Weyl/calibrated-coherent physical-volume
operator difference is less than 10^-200 on |u| <= 10^-2000.
The actual Weyl volume is greater than I/2 and within 10^-199 of I.

The corresponding full Hamiltonian operator ordering difference is less
than 10^1056. Exact two-propagator Duhamel comparison then bounds the
entire unitary difference in OPERATOR norm by 10^-944. No scalar phase,
implicit auxiliary term, generated spatial harmonic, matter/vector or
primitive contribution is discarded.

On the unchanged seed, the same-cutoff two-ordering readout mean differs by
less than 2 * 10^-200. The two explicit Weyl cutoff states differ by less
than 10^-943, and their volume means by less than 10^-199. The evolved
state's probability outside the coherent core is less than 10^-942, not
zero. The stronger S6.270 coherent-only cutoff bounds are not reassigned
unchanged to the Weyl comparison.

## Why this is an operator result

Direct integration of the defining Weyl kernel against normalized coherent
states fixes the cross-Gaussian prefactor pi^-d and the frame measure
(2pi)^-d dz. Product integration by parts in all 96 phase coordinates gives
the symbol-derivative weights [5,4,1] for orders [0,1,2]. Both full Schur
integrals give the dimension factor (pi/2)^48. The total weight is 10^96;
a safe explicit operator constant is 10^112.

All mixed derivatives with each coordinate order at most two, hence total
order through 192, are covered. The high phase-Cauchy argument uses radius
R/8 simultaneously in all 96 coordinates, inside the parent's complex
ball 4R. It does not upgrade the real-time profiles beyond their original
regularity.

The complete bump, transition, radial-partition and Leibniz estimates give

    J_n(A) = 2A 2056^n(n!)^3/R^n,  0 <= n <= 196.

Every one of the 196 actual growth ratios is less than 10^-9. With
D = Delta/4, the exact heat remainder is
-integral_0^1 t exp(tD)D^2 f dt. All derivatives required by Schur are
bounded by (96^2/32) J_4(A), for A = 10^1000 (interaction) or 10^-255
(volume error). Applying the explicit operator theorem, rather than the
symbol supremum alone, yields the quoted operator errors.

The Weyl-volume positivity is specific to this near-identity volume.
Independent positive-Gaussian-symbol counterexamples retain negative
Weyl excited-state eigenvalues; no generic positivity-preserving Weyl map
is asserted. The bounded real full Weyl Hamiltonian, its Schwartz kernel,
time continuity and exact Dyson evolution are treated in the written
proof. Both state and readout are conjugated in the physical picture.
The original translations and common zero-charge reducing sector remain.

## Source convention review

The author-hosted Didier Robert text,
[Propagation of Coherent States in Quantum Mechanics and Applications](https://www.math.sciences.univ-nantes.fr/~robert/proc_cimpa.pdf),
section 1 pp. 6-9, was consulted for the defining Weyl kernel and coherent
frame approach. The inspected printed Eq. (29) has a normalization factor
inconsistent with its diagonal unit-state check when combined with Eq. (28);
the later displayed frame inversion also omits its measure normalization.
Those prefactors are not imported. The exact Gaussian calculation and
independent quadratures derive the convention from the original kernel.

This review informed the independent normalization proof; the certificate
does not assume an unexamined external explicit operator constant.

## Validation evidence

The frozen packet has 18 source files and a 20-field native certificate.
Its exact audit contains 28 named checks, 32 scalar residuals, 43 proof
gates, nine controls, 235 rejected inputs, nine unchanged primitive rows
and 127 qualified matching rows.

Before freeze, Ruff fixed one unused import. Explicit changes then used a
TypeError for invalid amplitude types and bound the immediate quadrature
loop variables in their closures. All nine Python files pass Ruff check
and format check; all seven sibling-module export contracts resolve.
No scientific formula or ancestor was changed by these lint corrections.

- Private original-SymPy preflight: C0, 101.01734583405778 seconds, followed
  by all 872 independent scientific tests passing in 0.80 seconds.
- Accepted byte-identical repository/private/preflight copy: 18 files,
  90673 ASCII bytes, frozen 2026-09-14 20:51:22 UTC. All source hashes match,
  and every Markdown source has exactly one terminal newline.
- Fresh original-SymPy preflight: C0, 103.5385137089761 seconds, followed
  by all872 scientific tests passing in0.80 seconds. Every source hash
  matched the private preflight and accepted repository copy.
- Native full-ancestry report: C0 validation, 291117 ASCII bytes, 25 exact
  chunks and24 observed acknowledgments, protecting5533 scientific inputs.
  The raw un-reserialized SHA256 is
  970319d0dc460d46e9474e2ef9c17c1624aa7244b584461ac6cab1ff46771c7e.
  The private raw transfer was hashed before repository report creation;
  source18/field20/private/native validation passed. No chunk was truncated.
- Original-SymPy ordinary replay: The ordinary original-SymPy package replay completed with exit status zero: all 897 tests passed in 3179.88 seconds (52 minutes 59 seconds). Both captured test files were present and unchanged. Its zero exit was consumed at 2026-09-14 21:52:30 UTC.
- Original-SymPy CLI replay: The separate original-SymPy CLI replay completed with exit status zero, consumed at 2026-09-14 21:52:30 UTC.
- Captured full-P8 regression: the captured full-P8 regression completed with exit status zero: all
64887 tests passed in5339.17 seconds (1h28m59s). Its zero exit was
consumed at2026-09-14 22:57:17 UTC. The snapshot contained795 test
files with SHA256
cfbf372e8c3ed23ce7e2400271948aeee35f80442d1c2d2c9136b5ed2e6e07dc.
All captured files were present and unchanged. The runner exposed617
static namespace ancestors and5 guarded sibling helpers; its128
original-versus-adapter comparisons passed. Final exact-GCD counters
were68058 domain fallbacks,7388 exact descents and94 mixed fallbacks.
Only this full regression used the adapter; independent native,
private, fresh, ordinary and CLI runs used original SymPy.

The scientific suite covers independent defining-kernel and coherent
Gaussian quadratures, normalization and Schur integrals, all 193
product-weight coefficients, all 197 full radial partitions, all 196
growth ratios, positive-symbol Weyl counterexamples, Gaussian-radius
spectral comparisons, noncommuting unitary/readout tests and all rejection
controls. The expected ordinary suite adds 25 certificate tests for a
total of 897 cases.

The native, private, fresh, ordinary and CLI checks use original SymPy.
Only the separately audited captured FULL regression may enable the
exact-GCD adapter. Frozen ancestor files and reports remain unchanged.
The tests verify source-bound identities and finite inequalities; the
analytic operator arguments remain written proofs, not formally verified
proof-assistant theorems. Numerical fixtures are diagnostics and are not
the full P8 interacting evolution.

## Remaining original work

These comparisons concern the two finite regulated orderings only. They do
not establish the original interacting volume mean, an unregularized
singular Hamiltonian, a uniform mode/torus/cutoff-removal limit, physical
Wilsonian matching, omitted-loop bounds, all-energy UV consistency or
nonlinear global geodesic completeness.

The nine original primitive frontiers and prior matching rows are retained.
S6.261's physical binding remains REFUTED, the S6.265 integrability boundary
and S6.269 fixed-heavy-scale clarification remain, and scoped P8(a) is
unchanged. Original V/G/B and P8 are not closed.

The first full-regression attempt exited with status 4 during collection after 1247.91 seconds: the frozen S219 test could not import its existing sibling _full_spatial_reference module under isolated importlib collection. The complete-snapshot guard correctly rejected the missing test; this attempt supplied no scientific test verdict. No frozen source was edited. A harness-only guarded sibling-helper import context was added, with helper-byte hashes and origin/collision checks; all 59 focused harness tests passed, and the actual S219 helper imported from its exact source path. The retry started at 2026-09-14 21:28:01 UTC with the SAME 795-file snapshot and successfully collected all 64887 tests, with every captured file present and unchanged. The retry subsequently returned C0 with64887 tests passing; its final evidence is recorded above.
