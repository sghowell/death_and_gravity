# P8 S6.308: full-dimensional leading-ladder Coulomb phase

Native, original-SymPy ordinary/CLI and complete captured-suite acceptance
have passed. Original V/G/B/P8 remain OPEN. Scoped P8(a) is unchanged.

## Result and exact scope

This is a separately named leading straight-worldline ladder and
crossed-ladder class, with linear scalar stress and a free Einstein
graviton. It does not stand for the complete interacting amplitude.

On25/4<=s<=16,0<tau<=1 near either Bose endpoint, define
D=sqrt(s(s-4)),V=s^2-4s+2 and eta=V/(8pi*kappa*D), with the
unchanged kappa=1e800. In spacetime dimension4+2e, the numerator is
V_e=V+2e/(1+e). Its evanescent term must be retained during IR division.

After dividing by exp(i*eta/e) at each fixed perturbative order and
then taking e->0, the finite leading coefficients sum for|eta|<1 to
the Born pole times
C_eta(tau)=exp[i*eta*(ln(tau)-ell+2/V)]*Q(eta),
Q(eta)=exp(-2i*gamma_E*eta)*Gamma(1-i*eta)/Gamma(1+i*eta).

The general finite-difference argument cancels the fixed-order IR
poles. A separate Gaussian-Abel Fourier calculation recovers the same
finite phase. Its one-loop coefficient agrees with the full known
principal phase from S303, including the dimensional correction2/V.

For physical real eta,|Q|=|C_eta|=1. The bound|Q-1|<2eta^3,
together with eta<1/kappa on the whole stated energy window, bounds
the difference of the selected phase models by2e-2400 at every
positive transfer. Relative to the positive full Born, the comparison
costs at most a factor of two, giving4e-2400.

This estimate does not bound omitted recoil, nonlinear graviton,
hard-loop or matching contributions. The logarithmic phase remains
unexpanded. The limit order is fixed order, IR division, regulator
removal, and only then the finite generating-function sum; no
interchange with an unrestricted finite-regulator impact-parameter
integral is asserted.

Physical-edge unit modulus is not a complex Regge bound or a
forward holomorphic disk. The phase has nonzero complex monodromy.
No macroscopic classical condition Gm^2>>1 is assumed or satisfied.
The complete inclusive rate, actual quantum limit, common-parent
bounce and original V/G/B obligations remain open.

## Independent evidence

The general proof is calibrated by raw full-D rung formulas through
four exchanged rungs, finite IR-divided coefficients through six
loop orders, and the independent S303 one-loop phase. High-precision
Gamma evaluations check pole cancellation. Literal Gaussian-regulated
radial quadrature and a nonzero-transfer Abel/Riesz limit independently
test the Fourier transform. Separate controls detect monodromy and
the nonuniformity of a truncated logarithmic phase.

An initial scratch identity did not canonically combine positive-base
complex powers with principal-log exponentials. Making those positive
bases explicit fixed only the unfrozen representation; no branch
forcing, zero tolerance, backend replacement or frozen edit was used.

Cold70800 passed the whole original-SymPy preflight in90.578477s
and all503 science tests in2.74s. Independent fresh18747 forced
the complete S307 parent packets and passed in92.236052s, followed
by503 science tests in2.77s. All18 private/fresh/repository sources
were byte-identical before freeze.

The18 ASCII sources total66,486 bytes:308 named/scalar checks,
60 proof gates,8 controls,73 rejected inputs,9 primitive records,
164 matching records and6 historical qualifications.
The native helper received one specification, protecting6,221 inputs.
Its171,327-byte ASCII report arrived in15 chunks with14 ACKs.
SHA256:
b9698ae4dbd703253e7d6fdd011d00370fed6c72d256bfb246bc6a0ee216f1cf.
Both raw/source/count/frontier validators and the20-field REPORT AST
passed. Native raw freeze:2026-09-16 06:13:41 UTC.

## External acceptance

Original-SymPy ordinary72362 passed all528 tests in3297.90s (0:54:57).
Independent CLI86798 returned C0 with the scoped success verdict.
Both started2026-09-16 17:49:04 UTC. CLI completion was captured
2026-09-16 18:40:47 UTC and ordinary completion2026-09-16 18:44:19 UTC.
Complete outputs are retained.

Shared FULL40862 passed82,577 tests in 5644.89s (1:34:04), with all
869 captured files unchanged. Snapshot SHA256:
f0fcbc72c5385dca25337e1a08dbac522fc16f197223daa3208119cb55f18bfa.
The capture has691 namespace ancestors,5 helpers and128 adapter
contracts; counters are68,421 domain fallbacks,7,388 exact descents
and94 mixed fallbacks. It includes S308/S309, not S310 or later.
Completion was captured2026-09-16 16:25:14 UTC. The later FULL50300
also passed through S311:83,650 tests in5799.26s,873 unchanged files,
SHA256 bad08d33f93588b01f48977c40f9fd9ea2b9f8345f6966f1a4817653528931b9.

Native/direct/ordinary/CLI retain original SymPy; the audited exact-GCD
adapter is FULL-only. Publication checks the exact22-file scope,
frozen source/raw hashes and preservation of unrelated P4/P9 work.
Original V/G/B/P8 remain OPEN.
