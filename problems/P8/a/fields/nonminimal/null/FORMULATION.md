# P8-A.21 — a one-sided Wick-square cap does not rescue a null-line QEI

This is an explicit obstruction for an untransversely smeared null
segment in FOUR-dimensional Minkowski space. It is not a null
singularity theorem, a self-consistent SEE counterexample, or original
P8 closure. It keeps the physical nonminimal stress of A.20.

Use one real massless scalar, normal ordered against the Minkowski
vacuum, with 0<=xi<=1/4; xi=1/6 is the conformal member and xi=0 the
minimal known-answer member. Choose an affine sampling unit and
the null line gamma(s)=(s,0,0,s), |s|<=1, with ell=(1,0,0,1).
An arbitrary positive sampling length is recovered by dilation.

There is an explicit family of positive pure quasifree Hadamard
states omega_R, every member having finite total energy, such that

    <:Phi^2:>_R(gamma(s)) < -(3/5) A_R^2 <0,
    <:T_ab:>_R ell^a ell^b(gamma(s)) < -(3/10) A_R^2,
    A_R^2=hbar R^2 gamma0^2, gamma0>0, R>=1,               (1)

on the ENTIRE closed segment. At xi=1/6 the second coefficient
improves to 2/5. gamma0 is explicitly normalized from fixed
smooth compact wavepacket integrals in notes/state.md and
does not vary with R or the chosen sampler.

Consequently, for every fixed real nonzero g in C_c^infinity(-1,1),

    integral g(s)^2 <:T_ell,ell:>_R ds
        < -(3/10) hbar R^2 gamma0^2 ||g||_2^2 ->-infinity. (2)

Every state in this family satisfies the same one-sided
relative Wick-square upper bound w_R<=Phi_*^2 for ANY fixed
nonnegative Phi_*^2 on that entire segment, even Phi_*^2=0.
Thus no finite lower bound depending only on the fixed sampler
and such an upper cap can hold for this null average. In
particular A.20's one-sided state-amplitude hypothesis cannot
simply be transported from its timelike setting to this null setting.

The state is a squeeze of one ACTUAL normalized smooth compact
mass-shell wavepacket. Its fixed occupation n=1/3 and anomalous
covariance m=-2/3 satisfy n(n+1)=m^2, with positive quadrature
variances 1/6 and 3/2. The continuum mode, its commutator
normalization, all profile jets, and the full improvement
T_ell,ell=(partial_ell Phi)^2-xi partial_ell^2(Phi^2)
are retained. It is not a freely assigned pair of quadratic moments.

The family has no uniform total-energy, physical-momentum or
two-sided absolute Wick-square bound. Its total energy grows at
least as a positive constant times hbar R^2. At EACH finite R,
the COMPLETE-line averaged null energy is strictly positive:

    integral_R <:T_ell,ell:>_R ds
          =(4pi/3) A_R^2 I2 >0.                          (3)

This last statement is for this constructed family, not a new
general ANEC theorem. A compact negative average does not contradict it.
No limit state, cutoff-independent physical preparation, transverse
smearing bound, entropy bound or null SEE incompleteness claim follows.
