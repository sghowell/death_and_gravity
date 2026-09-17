# P8 S6.314: state transport and a two-real matched soft reference

Original V/G/B/P8 remain OPEN. Scoped P8(a) is unchanged.

## Result and strict scope

This successor supplies the state-transport term needed to connect
S309's one-real remainder dressing with S313's two-real subtraction.
Their counterterms are not identical and cannot simply be added.
The new expression reproduces the complete physical four-dimensional
zero/one/two-real tree densities with the original Born normalization,
recoil, parameters and polarization conventions.

Write P_sigma for the known leading-soft factor on the complete marked
state and P0 for its elastic counterpart. The required connector is

E1(x)=integral dB1(b)*[P_sigma_b(x-b)-P0(x-b)],

where dB1 is the nonnegative elastic one-soft Born measure. The
difference is taken before integrating db/b. Neither divergent
unsubtracted seed integral is assigned its own finite value.

The matched reference is P_match2=P0+D1+E1+D2. D1 is the unchanged
one-marked remainder dressing. D2 dresses the probability-level
two-real inclusion-exclusion remainder on its full two-marked state.
The Bose factor is retained: selecting two marked labels leaves
1/(2!*N!), and the seed measure already contains 1/2!.
The additional radiation shares the single remaining-energy cut
x-a-b, not two independent energy cuts.

For total marked energy R<=x<=1/8, the new estimates are

|Delta_sigma-Delta0|<5000*R*(1-ln R)/kappa,
TV(R2;x)<2e-725*x+2e-652*x^2,
|E1|/P0<1e5*x*(2-ln x)/kappa^2,
|D2|/P0<4e-725*x+4e-652*x^2.

Together with the unchanged D1 bound they give
|P_match2/P0-1|<1e-653, uniformly down to arbitrarily small positive
resolution x, with the relative difference tending to zero as x->0.
The resolution power is never expanded in a0*ln x. Pointwise
positivity of this named reference follows from P0>0 and the bound.

## Why the connector is finite and necessary

The two-real amplitude subtraction is B2=U+V-O. Squaring this is not
the same as the probability inclusion-exclusion subtraction:
|B2|^2-|U|^2-|V|^2+|O|^2=2Re[(U-O)conjugate(V-O)].
That retained interference is integrable. A second exact density
identity exposes the missing state-transport term when the one-real
and two-real probability remainders are combined. Its real term alone
still contains a single-soft logarithm, so its associated same-state
virtual difference is essential. E1 supplies that defined completion.

The finite angular conversion has a total-energy continuity modulus
R*(1-ln R), proved using nuclear norms and a joint radial Holder bound.
A uniform bound on the derivatives of the entire signed dimensional
soft sum then provides an integrable connector majorant at both energy
endpoints. The regulated contraction is not assumed positive; negative
regulated indices are explicitly retained. The regulator limit is
passed through the completed difference at fixed positive x, not
through separate divergent terms or fixed-real-multiplicity integrals.

Exact original recoil-current and full 434/47-tree configurations
calibrate the identities. Independent high-precision series and Beta
moment checks supplement, not replace, the general continuity and
matching proofs. The documented real-log integration substitution
resolves a scratch CAS branch result without discarding a physical
imaginary contribution or changing the backend.

## What is not established

This is a selected soft reference matched to stated real-tree sectors,
not a complete interacting detector distribution. Pointwise positivity
is not positive event-measure construction, monotonicity or unitarity.
No unknown finite hard loop or evanescent coefficient has been chosen.
Higher nonleading multiplicities, complete finite hard real-virtual
matching, the quantum state, absolute complex Regge control and the
common-parent bounce remain separate obligations. Original P8 is open.

## Immutable source and native evidence

Cold 14674 passed in 128.529s and passed 607 science tests in 0.20s.
Fresh 73211 forced the complete S313 parent packets, passed in 138.735s
and passed 607 science tests in 0.20s. All 18 private, fresh and repository
source files agreed byte for byte before freezing.

The 18 ASCII sources total 69,784 bytes, with 321 named exact residuals,
861 scalar entries, 131 proof gates, 8 controls and 73 rejected inputs.
All 9 primitive records, 170 matching records and 6 historical
qualifications remain unchanged except the explicitly appended item.

One native specification protected 6,335 inputs. Its 183,903-byte ASCII
report arrived in 16 chunks with 15 complete-chunk ACKs. SHA256:
ea01ed8699689c68db6d52ed9c9ebf25ee0015e478480e43391a60a68a03dd56.
Both raw/source/count/frontier validators and the 20-field REPORT AST
passed. Native raw freeze: 2026-09-16 18:34:41 UTC.

## External acceptance

Own ordinary 69538 passed all 632 tests in 3183.56s (0:53:03).
Independent CLI 98468 returned its explicit S6.314 success message.
Both used original SymPy and exited zero. Their completion was
captured at 2026-09-17 04:54:26 UTC, after the user's continuation.
Frozen source/native bytes remained unchanged.

The complete FULL 30261 passed 86,231 tests in 6121.67s (1:42:01),
with all 881 captured files unchanged. Snapshot SHA256:
c183015d4e2b956504d9ee5ed61c6880b7b8fdd671195d1cc68a5771b3411e8c.
It captured 703 namespace ancestors, 5 helpers and 128 adapter contracts.
Counters: 68,720/7,388/94. Completion was captured at
2026-09-16 21:26:12 UTC. It includes S314 through S315.

The newer complete FULL 55495 passed 87,609 tests in 6580.55s (1:49:40),
with all 885 captured files unchanged, 707 namespace ancestors,
5 helpers and 128 contracts. Snapshot SHA256:
9c307ad0e614657f3c68976825f9157af5060c2d86c6488f83417cb858b92924.
Counters: 68,723/7,388/94. Completion was captured at
2026-09-16 23:19:45 UTC. It includes S314 through S317.

The latest complete FULL 44112 passed 88,848 tests in 6791.09s (1:53:11),
with all 889 captured files unchanged. Snapshot SHA256:
6c99f6541e12a121814cad3c1fdf14a4f513b91c4c615c0efb0bd244364b7031.
It captured 711 namespace ancestors, 5 helpers and 128 contracts.
Counters: 68,749/7,388/94. Completion was captured at
2026-09-17 04:54:26 UTC. This snapshot includes S314 through S319.

Native/direct/ordinary/CLI use original SymPy; the exact-GCD adapter
is FULL-only. Publication validates the exact 22-file scope,
frozen source/raw hashes and staged bytes, preserving unrelated
P4/P9 work. Original V/G/B/P8 remain OPEN.
