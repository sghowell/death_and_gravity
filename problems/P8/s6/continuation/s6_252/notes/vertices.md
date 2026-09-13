# Whole coefficient vertices and independent response routes

Use the complete canonical H=-Omega A_c of S251, then substitute q=P^2/a^2 BEFORE differentiating its independent parameters(Jc,Theta,E,ell,A,Tcorr,a). The implementation records all seven first matrices and all28 symmetric second matrices. No A,Tcorr term, momentum cross term, normalization or potential is dropped. The physical-scale derivative includes -2q/a times the q derivative. Each TT Hessian is diag(a P^2,a^-3) in(h,a^3 h'), so its full variation is explicit too.

For the fixed state write W(t,s)=S(t,t*)[V0+iOmega/2]S(s,t*)^T. For two real symmetric full matrices A,B, the connected ordered Weyl-quadratic expectation is

C_AB(t,s)=Tr[A W(t,s) B W(t,s)^T]/2.

The final factor is an ordinary matrix transpose, NOT an adjoint and not a time-variable swap. The real part is symmetric noise. Twice its imaginary part equals the exact quadratic commutator susceptibility. In the initial canonical frame A_*=S(t,t*)^T A S(t,t*) and B_*=S(s,t*)^T B S(s,t*),

chi_AB=Tr[(A_* Omega B_*-B_* Omega A_*)V0]/2.

An independently varied impulsive Hamiltonian at s gives delta S(t,t*)=S(t,s) Omega B S(s,t*). Substitution into delta V produces the same chi. This proves the sign and full matrix order without a commuting-frequency approximation. The local second vertex gives the separate seagull. Normal ordering changes c-number means but does not erase that physical second variation.

For every finite real test combination the symmetric kernel is positive: its quadratic form is the expectation of the square of the centered Hermitian quadratic observable. The vacuum total oscillator energy has zero connected variance; replacing W^T by W^dagger fails this control. Full canonical basis changes transform vertices, state and propagators together and leave every scalar kernel invariant.

Independent coefficient probes still are NOT a physical metric/light variation. The latter needs the entire off-reference parent-to-coefficient map, density and canonical boundaries, nonlinear embedding contacts and gauge/constraint measure. Time derivatives in that map would act on the kernels and contact deltas with their proper adjoints. No such map is silently supplied by a collection of partial derivatives at the reference.
