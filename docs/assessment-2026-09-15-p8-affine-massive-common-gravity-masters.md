# P8 S6.283: common gravitational cut masters and external legs

All independent validation gates passed. Original V/G/B/P8 remain OPEN.

## Outcome

A single crossing-symmetric scalar-master basis reproduces the entire
original two-graviton normal cut and the entire pure-gravity part of the
massive elastic cut. Six ordered alternating-mass boxes appear once each,
although every box has both a massless and a massive cut. The full finite
evanescent trace and raw phase-space constants agree with S281.

The exact Feynman-parameter calculation supplies both triangle species
and both bubbles and proves the complete finite alternating box. The
literal off-shell scalar stress, harmonic graviton projector and scalar
self-energy give zero on-shell mass shift in D in the pure-gravity sector.
This is not a zero total mass shift with the retained heavy/matter loops.
The residue has distinct UV and IR pieces; four LSZ legs restore the
self term in S278's full massive analytic soft pole.

This is a cut-matched master representation, not the complete finite
rational amplitude or a physical-pole/local matching prescription.
Oepsilon bubble coefficients times UV poles, tadpoles and apparent Gram
denominators require further work. None is set to zero by assumption.

## Complete master dictionary

Use mu=m^2, Q=s-4mu, u=4mu-s-t and the exact S282 seed
M(a)=integral_0^1 dx/[4mu-a(1-x^2)], with the Feynman boundary.
The normalized scalar-loop convention is D=4-2eps with rGamma divided out,
rGamma=Gamma(1-eps)^2 Gamma(1+eps)/Gamma(1-2eps).
The original raw convention uses EP=-eps and retains
rGamma(-EP)(4pi)^(-EP)=1+EP(gamma_E-log4pi)+O(EP^2).

The complete alternating box obeys
I4(a,b)=4M(b)/a[1/eps+log(nu^2/(-a))]+O(eps).
Its literal Symanzik polynomial is
mu(alpha2+alpha4)^2-s alpha1 alpha3-t alpha2 alpha4.
The grouped-simplex and ratio substitutions reduce its full integrand to
r(1+r)^(2eps)/(A r^2+B)^(2+eps).
Exact radial and beta integrals give the displayed pole and finite term;
written uniform endpoint and tail bounds show the difference is O(eps),
not an omitted finite function.

The normalized massive triangle is
Gamma(1+eps)/(2eps*rGamma)*nu^(2eps)
times integral_0^1[mu-a*x(1-x)-i0]^(-1-eps)dx.
The finite two-massless-line triangle and both ordinary bubble parameter
integrals are explicit. All functions retain their i0 limiting prescription.

This independently derived box agrees with
[Ellis and Zanderighi, Box14](https://arxiv.org/html/0712.1851).
No four-dimensional massless approximation or potential/Born subtraction
is used as a substitute for the massive four-point amplitude.

Inside overall1/(16pi^2 kappa^2), each ordered(a,b) box has coefficient
V(b)^2, with V(b)=b^2-4mu*b+2mu^2. The three channel triangle/bubble
coefficients use the exact S282 polynomials A,B,R,D:
C00mu=2A(s)-6tu B(s)/Q^2 and B00=R(s)/30+tu D(s)/(10Q^2).
The massive triangle has coefficient-2V(s)H(s,t), where
H=-2s+6mu-2mu^2/s+tu/s.
Its massive bubble is the exact full regular polynomial moment.
VD=V+2mu^2 EP/(1+EP) and HD=H-2mu^2 EP/[s(1+EP)]
supply the box and massive-triangle finite evanescent coefficients.

Their complete normal cuts reproduce every S282 graviton term and every
S281 pure-gravity bare elastic finite term, not only their divergent parts.
A massless bubble-sum203/40*(s^2+t^2+u^2) is a calibration, not replacement
of the original external mass1.

## Literal external legs and remaining matching boundary

For the off-shell stress
T_mn=p_m r_n+p_n r_m-eta_mn(p.r-mu),
the Ward contraction is (p^2-mu)r-(r^2-mu)p.
Its full D-dimensional graviton contraction yields
2p^2 r^2+4mu p.r-2Dmu^2/(D-2).

After scaleless tadpoles vanish in DimReg,
Sigma(p^2)=raw[(4mu p^2-4mu^2/(D-2))B0(p^2;0,mu)-2mu A0(mu)]
divided by16pi^2 kappa.
The exact all-D Gamma relations give Sigma(mu)=0 in this sector.
Sigma'(mu)=mu[3/eps+7+3log(nu^2/mu)]/(16pi^2 kappa)
in the normalized convention. Its UV residue is4mu and its IR residue
is-mu; adding them prematurely would hide their different origins.

With inverse propagator p^2-mu+Sigma, four external legs multiply the
tree by1-2Sigma'. Their IR term supplies the missing massive self
contribution. Together with every crossed box and triangle pair pole,
this exactly matches S278's massive soft kernel including its self term.

This does not fix the finite physical residue, proper vertices, local
Wilson coefficients, tadpoles of the full four-point amplitude or any
higher-loop remainder. Other heavy/contact/matter sectors, detector
resolution errors, Regge input, the original quantum state/domain/measure
and bounce remain separate obligations. Scoped P8(a) is unchanged.

## Development and accepted source/native validation

The first assembled replay had all94 exact scientific residuals zero
but failed one structural-equality metadata gate; the unfrozen gate was
corrected to check the retained dimensional trace. Ruff then removed two
implicit module reexports, causing the first full preflight to fail on
import before scientific execution. Explicit inherited-module bindings
fixed this before acceptance. A temporary self-alias lint warning was
also resolved before the final source capture. No frozen input changed.

The final cold original-SymPy preflight passed in83.87922558397986s:
18 ASCII sources,66388 bytes,20 AST-derived report fields,7 export
contracts,94 named checks,97 scalar entries,36 gates,8 controls with73
actual rejections,9 primitive rows,139 matching rows and all6 historical
physical qualifications. All235 science tests passed in1.24s, exit0.
Exact formatted source bodies were captured and copied byte-for-byte.

The fresh repository preflight explicitly built all immediate-parent
packets before its own residuals. It passed in84.79416004102677s, and all
235 science tests passed in0.65s. All18 hashes matched the cold manifest.
Only then were the sources frozen.

The untouched original-SymPy native helper protected5746 unchanged inputs
before this packet and produced109968 ASCII bytes in10 complete chunks,
with9 acknowledgments. The last chunk starts at108000 and has1968 bytes.
SHA256:affd15ac3d06770b6c46b601138212bc1ecf537e91a636447141702568b72ab3.
Independent private and repository raw validators passed all byte/hash,
AST-field, residual/gate/control and retained-frontier comparisons.
The raw report is frozen. The helper returned READY without receiving
research probes or DONE;5747 inputs are protected including this report.

The original-SymPy ordinary replay exited0:260 tests passed in2878.64s
(0:47:58). Its separate CLI returned the complete common-gravity-master
and external-leg success verdict with original P8 OPEN, exit0.

## Completed independent replays and complete regression

The independent native, ordinary and CLI runs retain original SymPy GCD.
The exact adapter is restricted to the separately captured complete FULL
run. All captured ordinary files were present and unchanged.

The combined S282/S283 FULL exited0:71,564 tests passed in5493.59s
(1:31:33). All817 captured files were present and unchanged. Snapshot
SHA256:e282307b0fa171b2a4611e07636cdda3f6771869e68ff304d6683052ba325a2b.
The snapshot used639 audited static namespace ancestors and5 helpers;
all128 exact-adapter contract checks passed. Final counters were68130
domain fallbacks,7388 exact descents and94 mixed fallbacks.
This snapshot covers S282 and S283, not later S284/S285 work.
The rejected S279 archive remains outside the accepted test tree.

All independent acceptance gates passed. This certifies the stated exact
algebra and written proofs, not a Lean formalization, full interacting
physical matching or original V/G/B/P8 closure. The frozen source and raw
report bytes are unchanged; unrelated P4/P9 and scoped P8(a) are preserved.
