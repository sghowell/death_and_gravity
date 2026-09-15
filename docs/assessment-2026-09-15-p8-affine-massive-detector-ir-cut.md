# P8 S6.278: massive detector and original M1 cut

Source, private/fresh science, native, independent ordinary/CLI and complete
FULL regression validation passed. Original V/G/B/P8 remain OPEN; scoped
P8(a) and all earlier physical qualifications are unchanged.

## Outcome

The calculation supplies the actual original massless M1 contribution
to four-massive-scalar scattering at its first gravitational loop order,
together with a massive analytic soft kernel and an explicit subtraction
of all three crossed M1 low cuts. It is not a complete finite-gravity
amplitude, an exact quantum-vacuum theorem or a Regge remainder bound.

All three vacuum constants and the full fixed source functions remain.
The literal vacuum jets give canonical external tree mass squared1,
massless original M1, kappa10^800 and8piG=1/kappa. These are formal
coefficients around the retained flat quadratic reference, not exact
renormalized pole masses or a newly retuned quantum background.

Both stress-tensor Ward identities and an independent four-vector
contraction yield the complete production tree

    A(Phi Phi -> Psi Psi) = [s-(s-4m^2)x^2]/(4kappa).

Production has s+t+u=2m^2; the elastic four-Phi sum is4m^2. The distinction
is retained. Cartesian sphere moments and an independent spin0/spin2
sewing give, with the identical-intermediate1/2! and optical1/2,

    Im A_M1,s = [s^2-tu+2m^2s+6m^4]/(960pi kappa^2).

Optical positivity applies in the physical external region s>=4m^2.
The M1 threshold is0. Polynomial continuation below4m^2 is not a new
optical positivity claim.

The complete first M1 nonanalytic sector is
-c sum_z P_z log(-z/nu^2), c=1/(960pi^2kappa^2), modulo an unspecified
real crossing-symmetric local polynomial. Its forward transfer polynomial
is v^2+2m^4, so soft stripping does not remove the transfer logarithm.
The one-loop universal graviton-soft factor has no M1-species loop and
cannot cancel this species-specific contribution at that order.

## Massive soft kernel and the actual subtraction

The complete pair kernel is

    F(z)=[(z-2m^2)^2-2m^4]
         integral_0^1 dx/[4m^2-z(1-x^2)].

The integral proves analyticity away from[4m^2,infinity), including the
removable pseudothreshold0. Its exact coefficients and full four-leg/self
sum B=2(F(s)+F(t)+F(u))-m^2 retain the scalar mass. Conditional on full
soft factorization, changing physical detector resolution fromE0 toE1
multiplies the stripped amplitude by
exp[-B log(E1/E0)/(4pi^2kappa)]. Physical detector resolution is not the
dimensional regulator, and fixed-resolution decoupling differs from the
joint detector scaling limit.

The refreshed finite-apparatus framework of
[Bellazzini et al.,2512.13780v2](https://arxiv.org/html/2512.13780v2)
retains finite-coupling errors. Its numerical massless-scalar examples do
not furnish bounds for this massive spectrum. No numerical bound or
uncontrolled big-O remainder is imported.

For a stated subtraction boundary L>4m^2, the exact low piece is
-c sum_z P_z[log(-z/nu^2)-log((L-z)/nu^2)].
Keeping the real local polynomial untouched gives an analytic center
for this M1 sector. With D=L-2m^2 its forward v^2 coefficient is

    b20=-c[2log(D/nu^2)+log(L/nu^2)-12m^2/D-14m^4/D^2].

Its L derivative includes BOTH the dispersive s/u contribution and a
separate transfer term-c/L; dropping that term fails the exact identity.
At the original mass1, nu^2=1 and L=10^196, this nonlocal coefficient
has magnitude below2000/(960kappa^2), less than10^-990 of4lambda.
This is not a bound on the unknown local polynomial, other loops,
matching errors or a physical cutoff. The chosen L does not establish
Wilsonian validity up to that scale.

The finite remainder in the fixed-transfer route still requires actual
Regge data; see
[Tokuda, Aoki and Hirano](https://arxiv.org/html/2007.15009v2).
Light-mass-sensitive corrections in a different string-inspired model,
[Caron-Huot and Tokuda](https://arxiv.org/html/2406.07606v2), are not
transferred to this parent.

## Development and accepted validation

A preliminary private scratch probe first failed because its independent
pair-normalization comparison retained one extra beta. That scratch
comparison was corrected and passed; it is outside the scientific manifest.

The first assembled-source replay failed only the vacuum-constant
reference check: the independently specified old reference had omitted
evaluation atX=0. Bothu=0 andX=0 were supplied before source freezing.
The other18 source checks passed; the later three modules had not run.
The cut, soft and subtraction modules subsequently passed15,24 and14
named checks, with8,7 and7 gates, respectively.

The first whole private preflight passed in85.035330334s and all234
science tests passed in0.22s. Before acceptance, three multi-symbol
presence gates were strengthened to require every named symbol rather
than any one of them. Ruff fixed6 import blocks, formatted7 files and
left2 unchanged. No frozen source was edited.

The final formatted original-SymPy private preflight passed in
83.55202954099514s:18 ASCII source files,20 AST-derived report fields,
7 module export contracts,72 named checks,81 scalar entries,37 gates,
8 controls with72 actual rejections,9 primitive and134 matching rows.
All6 historical physical qualifications and all133 prior matching records
remain. All234 science tests passed in0.23s, exit0.

All18 accepted files total67985 ASCII bytes. Their exact formatted bodies
were captured, copied to the repository and compared byte-for-byte and
by SHA256 against the accepted private manifest. Sources were frozen at
2026-09-15 05:39:42 UTC and must never be changed in place.

The fresh repository original-SymPy preflight passed in84.72828233300243s,
then all234 science tests passed in0.35s. All18 source hashes equal the
private accepted hashes. Read-only Ruff checks passed for the9 private
and9 repository Python files.

The separate original-SymPy native helper protected5670 existing frozen
scientific inputs and rebuilt the complete ancestral report. Its84449-byte
ASCII result was transferred in8 exact chunks with7 acknowledgments,
without reserialization. The final chunk started at84000 and had449bytes.
Raw SHA256:
2cb6b7d65b559dcfc45dbaed1a8f058088cc15af87da32a9254ab70d9ca9afba.

Private-raw and repository-raw validations passed, checking all18 source
hashes,20 AST fields, all counts/gates/frontier records and all6 prior
qualifications. The raw report was accepted at2026-09-15 05:44:01 UTC.
The helper returned READY and was not sent any research probes or DONE.
There are5671 protected scientific inputs including this new report.

The independent original-SymPy ordinary replay actually collected259 tests
from both unchanged captured files and passed259 in2967.03s (0:49:27),
exit0. Complete output was captured without truncation; the final verdict
was consumed at2026-09-15 06:43:35 UTC.

The separate original-SymPy CLI replay returned exit0 with the complete
success line: original M1 gravitational cut and massive detector dictionary
replay passed, with original P8 explicitly OPEN. The verdict was consumed
at2026-09-15 06:40:55 UTC.

The shared S277+S278 complete regression started at2026-09-15 05:44:03 UTC.
It captured ACTUAL809 test files with SHA256
7516e55eff344c0107f309d3670d82c385c4cc06ffae276c7e1f3e4a9051b0fe,
631 static namespace ancestors,5 helper directories and128 exact
adapter-contract checks. Actual collection confirmed70401 tests and all809
captured files present and unchanged at2026-09-15 06:08:14 UTC.
ALL70401 tests passed in5482.56s (1:31:22), exit0; the complete final
verdict was captured at2026-09-15 07:15:48 UTC. Adapter counters were
68108 domain fallbacks,7388 exact descents and94 mixed fallbacks.
This snapshot excludes rejected S279 and later S280; it is not relabeled
as validating either packet.

## Remaining original task

The calculation fills actual missing IR data. It does not set the
two-graviton sector, other species, unknown local terms, higher loops,
finite detector errors or the Regge contour to zero. The original B
state/domain/measure and controlled common-parent matching remain open.
The finite classical endpoint ofS277 and corrected hybrid ofS275 retain
their separate, limited meanings.

Native/direct/ordinary/CLI use original SymPy. Only a captured FULL
regression may use the audited exact-GCD adapter. Publication requires
all actual gates, exact scoped staging and remote verification.
Unrelated P4/P9 work must remain untouched.
