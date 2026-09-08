# P8: first retained-vector local quantum coefficient

Original P8 remains open. The scoped photon objective, linear CD/M1
classification, physical matter frame and adopted V/G/B contract are
unchanged. This follows the [exact vector elimination/source remainder](assessment-2026-09-07-p8-exact-elimination-source-remainder.md).

## New result

[S6.47](../problems/P8/s6/matching/affine/kinetic/aligned/nonlinear/elimination/quantum/FORMULATION.md)
calculates the retained vector's first local one-loop coefficient at
fixed light fields. Its full four-component determinant keeps the
temporal constraint and all three massive-vector polarizations. An
independent physical longitudinal reduction agrees with the determinant.

The calculation states the Euclidean action, dimensional regulator,
continued transverse multiplicity, modified minimal-subtraction scheme
and renormalization scale explicitly. It derives the scalar Gamma
integral and finite vector constant, rather than treating a massive
vector as three four-dimensional scalar determinants. Independent
cutoff-logarithm and three-dimensional controls fix the pole sign.
The reference paper prints an opposite regulator definition, so no
sign-normalized match to its printed pole formula is asserted.

For a=1/gamma_t, b=1/gamma_s, m0^2=M^2/zeta_physical and
L=log(m0^2/mu^2), the finite local weight is

    B=2b^2[L+log(b)-1/2]+a^(3/2)b^(1/2)[L+log(a)-3/2].

Its pole weight is C=2b^2+a^(3/2)b^(1/2). On the clock, C=3 and
B=3L-5/2, but their actual first and second clock N derivatives do
not vanish. Thus source alignment does not remove the light-dependent
loop determinant. These are coefficient derivatives, not the full
quantum lapse Jacobian or stress tensor.

On the original closed coefficient tube, C<4 and, at mu=m0,
|B|<=3971/1296<4. The finite local potential is consequently smaller
than R^4/(144 Lp^2) times the reference scale M^2/tau^2, with
R=m0*tau and Lp=M*tau. R=1000, Lp=10^12 gives normalized
zeta=10^-6, m0/M=10^-9 and this particular local ratio below 10^-14.
That reference density is not a lower bound on a physical observable
which could vanish at the bounce. Finite counterterms and scheme/scale
dependence are not discarded.

## Remaining work

This is a constant-coefficient local calculation, not a global
Lorentz-invariant vacuum, a full curved-background quantum remainder,
a quantum-corrected bounce/spectrum or an interacting cutoff. V/G/B
remains open. No finite counterterm or frozen action was changed.

Work continues on the actual momentum-dependent vector mass-insertion
terms and a finite loop remainder after subtraction through fourth
derivative order. Curvature, mass-derivative, nonlocal/in-in and higher
loops, plus affine-complement/light/gravity contributions and matching
counterterms, remain separate. Small norm alone cannot protect the
original exactly saturated matter-cone condition; its sign/structure
still requires analysis. No user intervention is currently needed.

## Verification

The report pins 13 sources and fully rebuilds S6.46 and its ancestry.
It checks 18 named exact scalar identities, 16 continuous/interface
proof checks and 26 rejected inputs. The scientific suite passes
**21 tests in 3.00 seconds**; the ordinary suite passes **45 tests
in 226.37 seconds**, without the broad GCD adapter. The separate
seeded read-only CLI passes.

The full P8 regression through S6.47 passes **4599 tests in 830.38
seconds**, with no frozen checkpoint excluded. The unchanged exact
GCD runner passes 128 original tuple comparisons and records 5717
domain fallbacks and 6509 exact descents under the per-test seeded
recipe. The preliminary proof checker also caught a symbolic/native
Boolean interface mismatch, corrected before freezing without changing
the inequality. This is exact symbolic verification with written
scheme-explicit quantum and continuous proofs, not proof-assistant
formalization or external peer review.
