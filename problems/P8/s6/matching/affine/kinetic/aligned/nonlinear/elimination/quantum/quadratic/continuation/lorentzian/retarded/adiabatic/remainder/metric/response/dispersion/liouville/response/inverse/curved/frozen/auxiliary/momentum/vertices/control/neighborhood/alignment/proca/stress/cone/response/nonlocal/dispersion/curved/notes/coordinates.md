# Exact prepared time-coordinate channel

All stresses here use the Q normalization: rho and p denote
64*pi^2*L^2 times the fixed S6.82 dimensionless vector stresses.
Write S=rho+p. Thus rho'=-3H S. This rescaling changes no Ward identity.

For a prepared physical time change psi(u)=u+epsilon eta(u), identity
near the initial slice, pull back the minimal metric and ordinary Proca
field. On the homogeneous background N=psi', a=a0(psi). Each transverse
or constrained longitudinal canonical normalization obeys
g_N=g0(psi)/sqrt(psi'), while omega_N=psi'*omega0(psi). Consequently

    f_N=f0(psi)/sqrt(psi'), p_N=sqrt(psi')*p0(psi).

The prepared Cauchy data are exactly the original ones. There is no
independent change of state. For X=(|f|^2,Re(fp*),|p|^2), this gives

    delta X=eta X'+eta' diag(-1,0,1) X.

The code checks this directly in the actual new covariance evolution,
including delta d and delta omega. Substitution in all four normalized
physical readouts gives delta rho_Gamma=rho' eta,
delta p_Gamma=p' eta. The transformation is not merely a high-frequency
approximation.

All twelve new physical adiabatic source/readout variations obey the
same transformation order by order, with the varying frequency power
and a^-3 normalization retained. All six independently matched finite
local stress variations obey it as well; their three background Ward
identities pass. The source-dependent differences have the common
integrable limit proved in S6.85. This establishes the identity for the
actual renormalized prepared response, without appealing to symmetry
of a single-time in-in Hessian.

For general sources n=eta', v=w+H eta, linearity therefore gives

    delta rho_Gamma=R[w]+rho' eta,
    delta p_Gamma=P[w]+p' eta,

where R,P are the actual pure-scale response (n=0,v=w), not new response
prescriptions. Their varied Ward identity is

    R'+3H(R+P)+3S w'=0.

The already-fixed S6.82 scalar profile contributes delta rho_P=delta p_P
=S eta'. Its pure-scale normalized stress variation is zero. Retaining
these terms, transform the full background-cancelled currents
E_N=-delta rho_total, E_Z=3 delta p_total. Exact algebra gives

    Q_eta=S eta''+(S'+3H S)eta'-3H' S eta-3S w',
    Q_w=3P[w]+3p' eta+3S eta'.

Thus the time channel is local; no nonlocal lapse kernel is being
inverted or discarded. The first equation uses the varied Ward identity,
not an assumption that R[w] vanishes. The local eta/w cross terms are
weighted adjoints because S'+3H S=p'. Their compact density is

    -S eta'^2/2 -3H'S eta^2/2 -3S eta w'.

Both variations of this density are checked. The scalar 3P[w] still
includes the full selected-state exact-minus-adiabatic response and
its independently matched finite local operator. The fixed profile
cancels the coordinate-to-normalized current background contact; it is
not reset after choosing a source.

Finally the source map is a prepared bijection: eta is the zero-past
integral of n, and w=v-H eta. Its current transformation follows by
weighted integration by parts against a^3 du. With E_Z fixed, the
homogeneous kernel for E_N is E_N'+3H E_N=0. Prepared zero data force
E_N=0, so no original metric equation is lost.
