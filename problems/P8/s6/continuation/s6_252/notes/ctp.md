# The fixed-state full Gaussian closed-time-path functional

At a common finite regulator, choose independent REAL canonical modes, so a complex Fourier representation does not count P and-P twice. Use Weyl-ordered quadratic Hamiltonians H_+(t),H_-(t) and the fixed pure Cauchy covariance V0 from S251. The exact definition is Z=Tr[U_+ rho0 U_-^dagger]. The trace at the final time and the unchanged initial density are essential. This is not an in-out vacuum amplitude with a newly selected final state.

Choose a symplectic S0 with V0=S0 S0^T/2; the result does not depend on its orthosymplectic right-basis freedom. If S_+,S_- are the two full classical symplectic evolutions, put R=S0^-1 S_-^-1 S_+ S0. With its(q,p) blocks, define

alpha=(Rqq+Rpp+i(Rpq-Rqp))/2,
beta=(Rqq-Rpp+i(Rpq+Rqp))/2.

Their full CCR identities are alpha alpha^dagger-beta beta^dagger=I and alpha beta^T=beta alpha^T. In particular alpha is invertible. The Gaussian vacuum overlap is

Z=exp[-Tr Log(conjugate(alpha))/2],

where the logarithm is continued from identical histories ALONG THE ACTUAL HAMILTONIAN PATH. This follows by solving the transformed annihilation equations for the squeezed Gaussian and fixing its phase by the Schrödinger equation. Its modulus is det(alpha alpha^dagger)^(-1/4)<=1. A principal square root of an endpoint determinant is not a complete phase prescription: a one-oscillator2pi rotation has R=I but its metaplectic vacuum phase is-1. Equal CTP histories, in contrast, have R=I continuously and Z=1.

For completeness, let a=(q+ip)/sqrt(2), so U^dagger a U=alpha a+beta a^dagger. Inverting this canonical transformation gives U a U^dagger=alpha^dagger a-beta^T a^dagger. Consequently U|0>=c exp[a^dagger K a^dagger/2]|0>, with K=(alpha^dagger)^-1 beta^T=beta(conjugate(alpha))^-1. The inverse canonical identities make K symmetric and I-K^dagger K positive. The Gaussian norm fixes |c|=det(alpha alpha^dagger)^(-1/4), but not its phase.

Write a general Weyl Hamiltonian as a^dagger h a+Tr(h)/2+[a^dagger Delta a^dagger+a conjugate(Delta) a]/2, with h Hermitian and Delta symmetric. Its vacuum component gives c'/c=-i[Tr(h)+Tr(conjugate(Delta)K)]/2. Meanwhile alpha'=-i(h alpha+Delta conjugate(beta)), whence Tr[(conjugate(alpha))^-1 conjugate(alpha)']=i[Tr(h)+Tr(conjugate(Delta)K)]. The two equations prove the displayed continued logarithm starting from c=1. Thus the phase is obtained from the actual Weyl Schrödinger evolution, not inferred from the covariance. This derivation applies along the relative metaplectic path defining the CTP trace as well.

Keep any separate c-number action and its original formal grade separate. There is no new subtraction of its vacuum phase. The two independent TT oscillators use their complete two-by-two transformations. Their finite-regulator traces multiply the coupled scalar trace; the original H/Proca determinants are not removed or counted again as fundamental scalar/tensor loops.

For a Hamiltonian coefficient lambda_A, O_A=z^T H_A z/2 and Gamma=-i log Z. At equal histories,

Gamma_,A=-<O_A>,
Gamma_,AB=-<O_HAB> delta(t-s)-theta(t-s) chi_AB(t,s),

where chi=-i<[O_A(t),O_B(s)]>. The second coefficient derivative H_AB is the entire contact vertex, not an optional correction. With average/difference sources the quadratic imaginary part is+i/2 times the full symmetric connected noise. Direct overlap Taylor expansion through second order checks the mean, contact and noise with all mixing retained.

These formulas are time-retarded, not an assertion of spacelike microcausality for gauge-reduced variables. The finite Gaussian trace does not by itself remove a momentum cutoff or choose a physical curved counterfunctional.
