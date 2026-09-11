# P8: the free flat in/out state is Hadamard

Original P8 remains OPEN. This follows the
[complete free transition and out-energy bound](assessment-2026-09-10-p8-free-flat-dirac-in-out-energy.md).
No current research step requires user intervention.

## Same quadratic operator, now with a regularity theorem

[S6.166](../problems/P8/s6/continuation/s6_166/FORMULATION.md)
keeps the same free Minkowski Dirac operator, protected
SAT8 mass profile and in/out projectors. It establishes
Hadamard regularity by checking an identified external
theorem, not by treating finite particle energy as enough.

For s(x)=x/(1+x^8)^(1/8), an exact polynomial recurrence
controls every derivative. The complex-disk tail estimate
for x>=2 is

    |partial_x^n(s(x)-1)|<=5 n! 4^n x^(-8-n).

Oddness gives the other tail. In physical time,

    |partial_t^n(M(t)-M_asym)|
       <=5 Delta tau^8 n! 4^n |t|^(-8-n), |t|>=2tau.

These are all-order symbol bounds for fixed tau, not
bounds uniform in a sudden-transition tau->0 limit.

## Uniform one-particle in/out limits

The difference of the full and asymptotic Hamiltonians
is a bounded multiplication operator of norm |M-M_asym|.
The exact wave-operator integral gives, for T>=tau,

    ||W_in/out-W(T)||<=Delta tau^8/(56 T^7),
    ||c_in/out-c(T)||<=Delta tau^8/(28 T^7).

Both wave operators and their adjoints converge in norm;
their limits are unitary on one-particle L2 space. This
does not assert a global infinite-volume Fock implementer.

At the actual parameters and T=1, the bounds are below
1e-604 and 2e-604 respectively. These are operator-norm
errors, not local stress or Hadamard-seminorm errors.

## External theorem and convention dictionary

The external input is [Gerard and Stoskopf, Theorem 1.1](https://arxiv.org/pdf/2108.11955),
version 2 dated 18 October 2021. The checkpoint verifies
the relevant geometric, all-order mass-tail and spectral-gap
hypotheses and identifies the same limiting projectors.

Euclidean R^3 has bounded geometry; lapse one, shift zero
and spatial metric identity have vanishing asymptotic
differences. Both asymptotic masses are strictly positive,
and H_asym^2=p^2+M_asym^2 has a strict spectral gap.

The source's signature and evolution signs are not silently
identified with ours. Its Clifford dictionary is

    Gamma_0=-i gamma_+^0, Gamma_j=i gamma_+^j,
    beta_source=gamma_+^0, m_source=-M.

It yields exactly the physical Dirac operator. The source
Hamiltonian is minus the physical one, so spectral labels
exchange while the physical occupied projector is retained.
The theorem therefore applies to these same pure quasifree
in/out states.

This is an application of an external theorem, not a
repository proof or formalization of that theorem. It
does not supply quantitative local smooth-remainder
seminorms or the interacting curved state.

## Verification and immutability

The native report pins thirteen source/proof/test files
and twenty fields: 55 named identities, 55 scalar entries,
21 proof gates, 9 controls and 178 rejected inputs.

Final private science passed 356 tests in 21.05 seconds;
fresh repository science passed 356 in 20.89 seconds.
Ordinary replay passed 381 in 2119.64 seconds. Independent
native CLI replay passed. The complete captured P8 snapshot
passed 27620 tests in 3216.89 seconds, final exit code 0.

All 585 captured test files were present and unchanged.
Path-list SHA-256:

    ad07300f46bb1a0dac4e13b8a74105ffaa21d627d1de1e70ae0f7282e49fef75

The full-run exact GCD adapter passed 128 original tuple
comparisons. Final counters: 34090 domain fallbacks,
7340 exact descents, 88 mixed fallbacks. The snapshot
predates S6.167. Only full regression uses the adapter;
native, direct science, ordinary and CLI use unmodified
SymPy with interpreter-only allowances.

The 23221-character native report was transferred
losslessly in two chunks; all thirteen source hashes
were independently checked. Report SHA-256:

    b7843dda9b7b2a3687b1a5ce8ef2b57eb860f2517ad1a708e3115412a251df82

Private diagnostics initially found a structural matrix
comparison and an RK4 roundoff issue when tiny corrections
were repeatedly added to one. Before final verification
and freeze, the comparison was simplified exactly and
the integration variable changed to the deviation from
identity. The numerical tolerance was not relaxed.
No published scientific source or report changed.

Exact 17-file staging excludes later continuations and
all unrelated P4/P9 work.

## Remaining work

The sharper exact-frame transition/state-energy estimate
has passed its corrected private and repository science
and is awaiting the rebuilt native chain. A specified
absolute free one-loop local-energy calculation is in
private development, including its finite counterterms;
neither is counted as published here.

A free flat Hadamard state is not yet the interacting,
curved common-parent state. Absolute full stress, physical
cutoff and loop errors, quantum target matching, global
V contours/cuts and finite-gravity G remain open. Scoped
P8(a) and A.20-A.23 are unchanged; original P8(b), V/G/B
and original P8 remain OPEN.
