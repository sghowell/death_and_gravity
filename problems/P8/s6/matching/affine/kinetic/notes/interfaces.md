# Independent premises, physical units, and verdict boundary

The vector module reconstructs its four-vector Schur matrix from every
frozen quotient block and contracts the full 64-component stationary
solution. Its first variation differentiates the complete ADM metric,
inverse and Hessian before pulling back to the actual trajectory.
The scalar module separately contracts the old 64-component solution
and derives the coefficient/background-Hessian variation in a normal
frame. Their agreement is an executable interface, not an assumption
based on shared expected strings.

The bridges compare both reconstructions with the same functions
h=(1+u^2)^3, a=(1+u^2)^2, H=4u/(1+u^2), and with the original
J, Theta, Lambda, delta and matter w read from the pinned CD/M1
witness. In particular the scalar action receives the actual vector
mass matrix D_inverse=(8/3)diag(1,-1,-1,-1) on the trajectory and the
actual shifts alpha=3/(2h), d=21H/(8h), e=1/h.

The old complete small-symbol metric/free-matter Hamiltonian is rebuilt
from `p8.gamma.hamiltonian`. The scalar module's unmodified base action
is independently Legendre-transformed after the shift constraint and
its lapse is eliminated. The resulting Hamiltonians agree exactly,
including the matter momentum, matter/background mixing and lapse
source. This comparison is made before adding the vector, and does not
replace the new full scalar constraint calculation. The interface
using delta=(1+w/ell)/3 has ell!=0, as on the actual free-chi trajectory
at every finite time. All earlier source manifests remain frozen.

For the four-vector Schur action, the additional bound uses
47/100<p<11/20, proved from the closed target range
9/40<=p^2<=11/40. The polynomial p^3+2p^2+p-1 is increasing for p>0
and at the lower bracket equals 15623/10^6>0. Therefore

    D00 > 46869/1100000,
    -Dii > 14/55,
    ||D_inverse||_infinity < 1100000/46869 < 24.

These are continuous rational inequalities in the specified rest
frame, not grid observations. The looser p>9/20 bracket used safely
for a different bound in S6.37 would not prove this polynomial sign;
it is explicitly retained as a negative control in the vector tests.

## Nonunit normalization

Take physical coordinates x_phys=tau*xi and
V_normalized=tau*V_physical. The auxiliary mass term scales like
M^2*tau^2 while the curl-square action scales like its dimensionless
physical coefficient zeta_phys. After factoring the former positive
factor from the action, its coefficient is

    zeta_normalized = zeta_phys/(M^2*tau^2).

The isolated center vector block would have
m_phys^2=8M^2/(3zeta_phys) and
m_normalized^2=tau^2*m_phys^2. These are controls of the isolated block,
not the coupled spectrum or a certified heavy gap. For M^2=3, tau=2
and zeta_phys=5 the three values are 5/12, 8/5 and 32/5.
The independent unit implementations and direct action scaling agree.

The momentum variable is q=(tau*k_physical)^2. The limiting negative
scalar-pivot threshold as u approaches zero is q=3597/3200; in the
nonunit example this is k_physical^2=3597/12800. The physical kinetic
proof uses u!=0, not an evaluation of the singular unitary chart at
the center. None of these numbers specifies an eigenfrequency or a
low-energy EFT cutoff. In particular zeta_normalized -> 0 removes the
extra vector pair exactly at zero and cannot be used as a regular
three-mode kinetic-rank limit.

## Frozen-parent replay

The wrapper checks the S6.37 report hash and fully rebuilds and validates
that parent report, including its S5.5 parent replay. That process also
checks the older source manifests and original formulation/contract
pins; it does not recursively rerun every older scientific report.
The complete P8 regression is a separate verification step.
New source hashes cover all local code, tests, formulations and notes.

Original P8 remains open. HOMOTHETIC is a spectator result, while
TRACE_FORM and QUOTIENT_VECTOR have their own precisely scoped health
rejections. No conditional gravitational-dispersion assumption has
been replaced by a flat-space positivity formula, and no missing UV
matching estimate has been declared solved by an auxiliary inverse.
