# S6.72: causal inverse of the exact isolated massive dispersion block

For fixed h>0 and m>0, retain the two physical lapse/log-scale
currents, the S6.67 finite fourth-order local coefficient F, and
the exact S6.69 thrice-subtracted *flat vacuum* bubble. Divide
their sum by s^4/(64*pi^2). The resulting matrix B(p), p=s^2, is

B(p)=F-p/4 integral_0^1 y^2 M(y^2)/(4m^2+(1-y^2)p) dy.

No leading-log replacement, tree term, zero/second-order local
term, curved mode kernel or new finite subtraction is included.
The physical choice h=(1+u^2)^3 lies in [1,125/64] on I.

## Scoped result

B has no determinant zero on its first p sheet, including its
finite threshold and both open cut banks. Its inverse R has the
convergent matrix representation

R(p)=Rinf-integral_(4m^2)^infinity rho(tau)/(p+tau) dtau,

with rho strictly positive definite on the open cut. Rinf is a
nonzero rank-one matrix. The corresponding causal operator is
Rinf times the identity plus convolution by a real matrix
Kreg in L1(0,infinity). In particular it is a bounded C0
endomorphism on every finite interval and on bounded continuous
half-line inputs. On zero-past distributions it is a two-sided
inverse of this isolated B, not of the full response.

The regular kernel has an integrable logarithmic small-time
singularity and an integrable threshold late-time tail. The
frequency sine integral is oscillatory, not absolutely convergent.
An exact positive static spectral moment bounds the primitive
and supplies an explicit prepared C1-to-C0 estimate. No numerical
L1 norm of Kreg is claimed.

A smoothly varying h can be put inside and outside the frozen
convolution to define a bounded causal preconditioner. This
definition is not an identity for the actual curved response.

## Evidence and boundary

The report fully rebuilds S6.71, which rebuilds the preceding
physical response and dispersion checkpoints. Algebraic endpoint,
sign, cut, Gram and moment checks support the self-contained
written analytic proof in notes/first-sheet.md,
notes/causal-kernel.md and notes/next-curved.md.

The pre-existing selected state and physical finite prescription
are unchanged. This flat-vacuum inverse is not their full curved
response inverse. Spatial response, full tree-plus-loop feedback,
quantum stability/cones, interactions, cutoff and common-parent
V/G/B remain open. Original P8 is not closed.
