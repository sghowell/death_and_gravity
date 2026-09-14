# Complete finite-band generator difference and energy estimate

Start with the entire S254 current coupled Hamiltonian, retaining
the full heavy mass multiplier, both heavy charges, heavy lapse
source, temporal-vector contact, all lower Euler/Ward contacts,
and the original -H_hat Qb Pb term. The current unforced equations
set the lower Euler contacts to zero only after the full Hessian
has been derived. The fixed reference has L0=3 Tcorr and Vvv=9 A/2;
these tiny but nonzero terms remain in the comparison.

Let R_ref be the fixed symplectic shear from notes/reference.md.
On the coupled eight phases use the CONSTANT mass-adapted weight

    W=diag(P,P,sqrt(P^2+mH^2),P/sqrt(zeta),
           1,1,1,sqrt(zeta)).

Each tensor and transverse-Proca pair uses diag(P,1).
The complete16-phase transformed generator is

    A = W R_ref^-1 Omega H_current R_ref W^-1
        -W R_ref^-1 R_ref_dot W^-1,

with identity shear blocks on the remaining pairs. The second term
is nonzero and is retained in each evolution equation. It cancels
only in the difference against the same fixed-reference equation.

Every entry of every partial derivative with respect to all24
independent current coefficients is bounded by exact outward
rational intervals over the entire band[10^64,2*10^64] and the
displayed convex coefficient box. The heavy mass lies in
[10^197,10^198]; its energy weight is kept in[10^98,10^100].
The all-entry sum, which bounds the Euclidean matrix operator norm,
is less than10^180. This is an explicit full finite-band bound, not
a principal-symbol estimate or a large-momentum extrapolation.

The trajectory proof puts every current and fixed-reference
coefficient, and every straight segment between them, inside this
box. All24 coefficient differences are smaller than10^-143.
The reference connection has zero derivatives with respect to
those current coefficients; all6144 such scalar derivatives are
checked. It follows that

    ||Delta A|| <10^180*10^-143 =10^37.

The reference energy E=x^T Q_ref(t) x/2 has eigenvalues between
10^-12 and10^12. Its extra root-energy rate under Delta A is
bounded by sqrt(cond Q_ref)||Delta A||<=10^12||Delta A||.
Together with the original reference energy identity this gives

    |(sqrt E)'| <= (10^28+10^49) sqrt E.

This holds almost everywhere, with the corresponding Dini
interpretation at zero. On the FULL slab length2*10^-60, the
exponential is bounded by its positive-series majorant

    exp(2 T rate) <= 1/(1-2 T rate) <2.

The estimate works in either time direction and between any two
times in the slab. It applies separately to real and imaginary
Fourier parts and hence to the full-band L2 reference energy.
There is no time integration of the actual enormous oscillator
frequency in the executable proof.

The complete original current canonical map is retained at both
endpoints, including the longitudinal scaling and central
symplectic map. Exact inverses are checked. Direct entry-sum
interval bounds give each conversion norm<10^150.
The reference metric adds at most10^12 to the Euclidean norm
conversion, so the full original unweighted phase propagator is
bounded by2*10^312<10^313. That large ceiling is intentionally
reported alongside the small reference-energy growth.

The existing S256 epsilon=10^-6 finite-growth example concerns a
different, much larger displacement. It is not refuted by this
epsilon<=10^-230 family. Neither estimate establishes a physical
Wilsonian cutoff, all-momentum control, a nonlinear inhomogeneous
Cauchy theorem, or a self-consistent interacting quantum mean.
