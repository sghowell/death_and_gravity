# Convergent regulated dynamics with the original free flow

Let a(u)=(1-D)g_ext(u), with the entire local symbol, moving-chart
connection, reference-flow subtraction and cutoff contacts
specified in the source and ordering notes. On the common
compact time interval a is a scalar plus a compact smooth real
function, norm-continuous with all needed derivatives. Put

    A(u)=Q_V(a(u)),   B(u)=||a(u)||infinity.

Then A is bounded self-adjoint and norm-continuous, and
||A(u)||<=B(u)<=||g_ext(u)||infinity+||D g_ext(u)||infinity.
All statements initially concern the whole L2 Hilbert space.

Iterate the integral equation U_I(t)=I-i integral_0^t A(s)U_I(s)ds.
If Lambda(t)=integral_0^|t| B(s)ds, the ordered k-simplex integral
has norm at most Lambda^k/k!. Thus the Dyson series converges in
operator norm, uniformly on compact time intervals, and solves
the equation. The same construction backward from an endpoint
gives its inverse. Differentiating U_I^*U_I and using A=A^*
gives exact unitarity and the composition law.

An improved truncation estimate follows only AFTER unitarity
has been established. Iterate Duhamel M+1 times but leave the
exact endpoint propagator in the final ordered integral.
Its norm is1, giving

    ||U_I-sum_{k=0}^M D_k|| <= Lambda^(M+1)/(M+1)!.

No extra exp(Lambda) is needed. This is an operator norm
estimate. It must not be confused with the separate sup-norm
Weyl-symbol calibration estimate, which is not on its own a
Weyl-operator norm bound.

Let U_ref be the metaplectic propagator of the full original
quadratic reference h_ref. The actual regulated propagator is

    U_full=U_ref U_I,
    H_reg=H_ref+U_ref A U_ref^*.

Quadratic real Hamiltonians on finitely many canonical pairs
have their self-adjoint metaplectic generators with Schwartz
core. This can also be seen by the dense finite Hermite span
of analytic vectors: repeated quadratic operations change
total occupation by at most2 and their factorial bounds give
a positive analytic radius for every finite Hermite vector.
For smooth time-dependent quadratic coefficients the smooth
symplectic flow has its continuous metaplectic lift.

At each time the second term in H_reg is a bounded self-adjoint
perturbation, so H_reg is self-adjoint on the reference domain
and retains its Schwartz core. More is needed for a strong
time-dependent Schrödinger equation, and it is available here.
The compact smooth part of the coherent symbol has a kernel
which maps L2 to Schwartz: every configuration derivative and
polynomial weight of the Gaussian displaced window is uniformly
bounded over its compact phase support. Cauchy-Schwarz bounds
the bra by the L2 norm. These give a finite L2-to-each-Schwartz-
seminorm operator bound, uniform on compact time intervals.

Separate the retained scalar background phase. In each remaining
Dyson term use one outer compact-part factor for the Schwartz
seminorm and the remaining factors for the operator norm.
The resulting factorial majorant proves convergence in every
Schwartz seminorm. Thus U_I and U_full preserve Schwartz and
the displayed full Schrödinger equation holds there strongly.
The large heavy frequency is not removed; finite mode number
and finite time, rather than a small-frequency approximation,
are what make this construction available.

The original compact three-translation group preserves the
seed, the full reduced source, the free flow and the chosen
cutoff. Coherent covariance gives T Q_V(a)T^*=Q_V(a composed
with the inverse phase action)=Q_V(a). Therefore A,U_I,U_ref
and U_full commute with all group unitaries. The common
zero-charge subspace is reducing and contains the SAME psi0
already, without conditioning or a new normalization. The
nonzero local constraints were solved classically; this is
not a full continuum diffeomorphism or BRST regulator.

For the positive coherent volume F_ext define

    O_reg(u)=U_ref(u) Q_V(F_ext(u)) U_ref(u)^*.

The positive finite regulated readout is exactly
<U_I psi0,Q_V(F_ext)U_I psi0>. If U_M approximates U_I in
operator norm with error epsilon and B_F>=||F_ext||infinity,
expanding the two sesquilinear forms gives the error bound

    |exact mean - approximate mean| <= B_F(2epsilon+epsilon^2).

A calibrated volume uses a distinct operator and its separate
conditional positivity criterion. No original unregularized
interacting mean, evaluated P8 value, quantum matching cutoff
or regulator-independent limit is inferred.

Independent diagnostics use a TWO-configuration-oscillator
coherent integral, compressed only for testing to total
occupation<=3, dimension10. The two real symbols are

    [(q^2-p^2)+(q^2+p^2)^2/10] exp(-|alpha|^2/2),
    [(q dot p)+(q dot p)^2/20] exp(-2|alpha|^2/3).

They generate Hermitian noncommuting matrices H0,H1 and both
commute with J=i(a1^*a2-a2^*a1). The positive test symbol
[(q^2-p^2)^2+1]exp(-|alpha|^2/2) gives a strictly positive
compressed matrix. The exact coherent monomial formula
integrates over the ENTIRE phase space; it is not a finite
CCR construction. A total-number cutoff preserves rotations,
whereas no finite matrix can have a commutator of trace equal
to the nonzero trace of iI.

For H(t)=H0+tH1, exact polynomial Dyson coefficients retain
time ordering and satisfy the differential equation and
unitarity identities order by order. An independent numerical
matrix ODE comparison checks the factorial bound, unitarity,
translation charge and readout error. This Gaussian-decay
fixture is not the actual compact P8 cutoff, full P8 symbol,
nonlinear evolution or a physical error estimate.
