# Actual ordinary-Proca addition at zero vector background

The S6.81 physical vector action is

    sqrt(-g) [-zeta F(W)_{mu nu} F(W)^{mu nu}/4 - W_mu W^mu/2],
    zeta=10^-6.

The CONSTANT rescaling A=sqrt(zeta) W gives the canonical Maxwell term
and m^2=1/zeta=10^6. The native identities keep this factor and agree
with the actual parent mass. There is no time-dependent field rescaling
or dropped derivative interaction in this step.

At W=0 the vector action and its first vector variation vanish for every
metric. With vector perturbation epsilon A and metric perturbation eta h,
it starts at epsilon^2; metric-vector mixing starts at eta*epsilon^2.
proca.py varies all ten independent inverse-metric first jets, including
the first volume variation, and checks the mixed eta/epsilon Hessian is
zero. In the local physical inertial frame (-,+,+,+), the quadratic form is

    (E^2-B^2)/2 + m^2 A_0^2/2 - m^2 A_i^2/2.

The scalar/metric quadratic Hessian is therefore unchanged, including
the lapse and shift rows. Linear gravitational constraints receive no
vector source: the vector stress starts at quadratic order. Eliminating
the already regular light constraints cannot introduce a quadratic
light/vector cross term. The vector temporal constraint acts inside its
own block. This argument is at the actual zero-vector classical solution;
it does not transfer a quantum stress tensor or state to a new solution.

## Three physical modes, not four unconstrained components

Freeze coefficients in a local physical tangent frame to read the
principal characteristics. For one real longitudinal Fourier quadrature
the action is

    L_L=(v_dot-k A_0)^2/2 + m^2 A_0^2/2 - m^2 v^2/2.

The temporal equation gives A_0=k v_dot/(k^2+m^2). Substitution yields

    L_L= m^2 v_dot^2/[2(k^2+m^2)] - m^2 v^2/2.

The finite-momentum kinetic coefficient is positive for finite m>0,
and omega^2=k^2+m^2, as for the two transverse modes. The front speed is
one at fixed finite mass. Equivalently the actual covariant massive
vector equation and its divergence constraint give three physical wave
polarizations on the physical metric cone. The raw four-component
unconstrained principal determinant is singular and is not counted as
four healthy modes. A time-independent FLRW dispersion or uniform
infinite-mass limit is not asserted.

## Finite spectator family

For n independent ordinary Proca species, with zero classical vector
backgrounds and arbitrary fixed finite positive masses, the physical
principal factor, up to nonzero normalization, is

    (c^2-s^2)(1-s^2)^(3+3n).

There is one clock mode, one matter mode, two tensor modes and three
physical modes per vector. The native 2-by-2 mass-mixing anchor also
checks that finite lower-order masses do not change a second-order
principal coefficient. The all-n product follows from the decoupled
blocks, not from extrapolating three sampled species counts.

n=1 is the actual existing ordinary-Proca candidate; n counts TOTAL
species of this separately specified classical family. n=0 is its
vector-free classical member, not a statement that the old quantum
Proca state or background expectation values can be removed unchanged.
Adding any finite n leaves the original faster scalar factor. Pure
classical spectator addition therefore does not repair this cone.
