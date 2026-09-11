# Out-particle energy and local in/out energy difference

Keep the symmetric Dirac stress tensor, treating M as a
prescribed external scalar under metric variation. Its
00 component is i psi^dagger bidirectional_partial_t psi/2.
The equations of motion reduce its homogeneous expectation
to the one-particle Hamiltonian bilinear.

For context, see
[del Rio et al., section IV, equations40-41](https://arxiv.org/pdf/1703.00908)
and the discussion following42 of state-independent UV terms.
All constants and state-difference bounds below are derived here.

Set C=2^42 and E=sqrt(p^2+m0^2). The complete new transition
bound is |beta|<=C pDelta/(tau^4 E^6).
The asymptotic out-particle energy from the in-state remains

    rho_particles=4N/(2pi^2) integral p^2 omega_out |beta|^2 dp,

with N=6 active color/flavor components, two helicities,
and particles plus antiparticles. Inert constant-mass
flavors have no transition. Using omega_out<=2E and

    integral p^4/E^11 dp=8/(315m0^6)

gives the sharper bound

    rho_particles<=32NC^2 Delta^2/(315pi^2 tau^8 m0^6).

At the actual parameters, using pi^2>9, this is approximately
1.25227336376e19, strictly below 1e20. It strengthens rather
than changes the older looser certificate.

Now compare the two Hadamard states on the SAME free
time-dependent operator. Write Delta P for the difference
of their equal-time occupied covariances. For one helicity,
two pure rank-one projectors differ with eigenvalues
+/-|beta|. Both helicities give trace norm4|beta| per
color/flavor. Evolution by the same Hamiltonian preserves
this trace norm at every real time.

With one identical local subtraction and finite counterterm
prescription, state-independent terms cancel in this
difference. The Hadamard property supplies the proper
local definition; the new bound also makes the relevant
momentum integral absolutely convergent. Since
||H(t)||=sqrt(p^2+M(t)^2)<=2E,

    |rho_in(t)-rho_out_state(t)|
       <=4N/pi^2 integral p^2 E |beta| dp.

The radial integral

    integral p^3/E^5 dp=2/(3m0)

therefore proves, uniformly in real time,

    |rho_in(t)-rho_out_state(t)|
       <=8NC Delta/(3pi^2 tau^4 m0).

The actual rational allowance is approximately
2.36931798578e410, strictly below 1e411. Its ratio to
the named kappa*m_Phi^4 reference is below 1e-389.
The out-particle allowance's reference ratio is below
1e-780. Neither is a relative error against bounce
density, which is zero at H=0.

The local comparison retains covariance coherence terms
linear in beta. It does not replace them by occupation
numbers |beta|^2. Nor does it set the out-state's local
energy to zero during the mass transition.

The fermion subsystem exchanges energy with its external
mass source. At the state-difference level,

    partial_t Delta rho=M'(t) Delta<bar psi psi>.

This follows by differentiating tr(H Delta P): the
commutator evolution of Delta P has zero trace against H,
leaving only H'=beta_D M'. The bilinear difference integral
is also convergent under the p^-5 transition bound.
No isolated conservation law for the time-dependent
fermion subsystem is imposed.

An absolute stress estimate would additionally need
the local vacuum-polarization and finite counterterm
contributions. A curved/interacting parent requires still
more. Neither is supplied by a difference of free states.
