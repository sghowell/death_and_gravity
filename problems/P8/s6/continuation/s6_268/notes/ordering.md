# Normalized coherent quantization and full ordering contacts

Use unit CCR after the original physical normalization. Let psi0
be the SAME full pure finite Gaussian reference with covariance V.
For z=(q,p), let D(z) be Weyl displacement and psi_z=D(z)psi0.
Define the weak bounded-symbol integral

    Q_V(a)=(2pi)^(-d) integral a(z)|psi_z><psi_z| dz.

Its normalization can be proved without a symbol convention:
in configuration space integrate the p-dependent exponential in
psi_z(x)conjugate(psi_z(y)) to get delta(x-y), then integrate
|psi0(x-q)|^2 over q to get1. Equivalently the coherent analysis
W psi(z)=(2pi)^(-d/2)<psi_z,psi> is an isometry into phase-space
L2 and Q_V(a)=W* M_a W. Consequently Q_V is positive and unital,
real bounded a give bounded self-adjoint operators, and

    ||Q_V(a)|| <= ||a||infinity,
    Q_V(a)*Q_V(a) <= Q_V(|a|^2).

The second inequality follows by inserting the orthogonal
projection WW* between the two multiplication operators.
These are infinite-Hilbert-space statements; no finite matrix
satisfies the exact canonical commutation relations.

The Weyl symbol of a rank-one projector is (2pi)^d times its
normalized Wigner density. Integrating the displaced rank-one
symbols therefore gives Q_V(a)=OpW(a*w_V). For the Gaussian
reference w_V has its ENTIRE covariance V, so with

    D=(1/2) sum_AB V_AB partial_A partial_B

the Weyl symbol is e^D a. The same-state expectation adds a
second Wigner convolution:

    <psi0,Q_V(a)psi0> = (a*w_(2V))(0) = e^(2D)a(0).

This is the Husimi density of psi0 relative to its own displaced
windows. Its covariance is2V, notV. Every q-p, field-field and
tensor/matter covariance in the full original preparation remains.

The primary convention comparison is
[de Gosson, section3](https://arxiv.org/pdf/1907.02471):
Proposition3 uses unnormalized Lebesgue measure; dividing by
(2pi hbar)^d gives the convention above. Only that convolution
and Corollary4's window covariance are used. The complete
relevant PDF pages6-8 were inspected; no later mixed-state
theorem is imported. All estimates here have the direct
normalized proofs just given.

For the interaction choose the explicitly named first-Weyl
calibration a_match=(1-D)g_ext. If g_ext is a constant plus a
compact smooth symbol, every derivative below is bounded.
Differentiate e^(tD)(1-tD)g_ext to obtain

    d/dt [e^(tD)(1-tD)g_ext] = -t e^(tD)D^2 g_ext.

The Gaussian heat semigroup is a sup-norm contraction. Integration
from0 to1 gives the exact identity and estimate

    e^D(1-D)g_ext-g_ext = -integral_0^1 t e^(tD)D^2 g_ext dt,
    ||Weyl-symbol error||infinity <= ||D^2 g_ext||infinity/2.

This is a SYMBOL bound, not automatically a Weyl-operator norm
bound. Positivity of a Weyl symbol likewise does not by itself
prove positivity of its Weyl operator. The bounded coherent
operator itself has the separate bound
||Q_V(a_match)||<=||g_ext||infinity+||D g_ext||infinity.

The calibration is exact on cubic polynomials because D^2 then
vanishes. On quartics the remaining symbol is -D^2 g/2. The
full correlated/squeezed four-phase-coordinate fixture gives
D^2 g=826341/31360, explicitly nonzero. Thus this prescription
is not exact nonlinear Weyl ordering or an implicitly justified
counterterm of the original theory.

All derivatives of the cutoff are mandatory:

    D(chi f)=chi Df+f Dchi+sum_AB V_AB (partial_A chi)(partial_B f).

Using chi(f-Df) in place of (1-D)(chi f) discards both last
contacts and defines a different operator. Comparing two
calibrated cutoffs must include those derivative differences.

For the physical observable retain the complete implicit lapse
root and nonlinear canonical-to-invariant map:

    N_A=N_i z^i_A,
    N_AB=N_ij z^i_A z^j_B+N_i z^i_AB.

All twelve first derivatives and all144 second invariant
derivatives are inherited from S266. For F=exp(3v)U(N),

    F_AB=exp(3v)[
      9U v_A v_B+3U_N(v_A N_B+v_B N_A)+U_NN N_A N_B
      +3U v_AB+U_N N_AB].

Both nonlinear second contacts, the mixed contractions, and the
actual U=R_full^(-3/4) remain in D F. A calculation using only
U_NN Var(linear n) is insufficient.

The uncalibrated coherent volume Q_V(F_ext) is positive.
A calibrated volume Q_V((1-D)F_ext) has this positivity proof
only if ||D F_ext||infinity<inf F_ext. That quantitative gate
is not evaluated here. These two readouts and their ordering
difference are kept distinct; neither is promoted to the
original interacting physical-volume expectation.
