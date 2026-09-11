# Geometry, units and the physical quadratic operator

The earlier analytic affine-domain target inherits
p8_exceptional_vacuum.family.original_retuned_tree_scalar.
That tree reads p8_affine.dictionary.target(), which reads the pinned
witness-CD_matter.json. Its scale factor is (1+t^2)^2, not the D-only
scale factor 1+t^2. The added retuning is proportional to (X-1)^2;
its background value and first X derivative vanish.

Use the existing mass-one reference time unit and geometric bounce
scale one. The canonical clock is Phi=sqrt(kappa)t. R=1e300 and
sqrt(kappa)=1e400 give tau=1e-100 for the SAT8 mass transition.
Tau is not the geometry's bounce time scale.

There are six active color/flavor copies with masses
M=m +/- Delta s(t/tau), s(x)=x/(1+x^8)^(1/8), and 36 constant-mass
copies with M=m. Here m=1e200, 0<=Delta<=3e197, m0=.99m.
The same estimates cover every copy conservatively using Delta's upper
bound. Setting Delta=0 removes the varying-mass contribution but not
geometric mixing.

In signature(+---), the physical curved equation is

    [i gamma^0(partial_t+3H/2)+i gamma^j partial_j/a-M]psi=0.

The canonical Cauchy variable psi_c=a^(3/2)psi has the flat spatial
L2 norm. Its helicity Hamiltonian is q sigma1+M sigma3, q=p/a,
omega=sqrt(q^2+M^2). The negative asymptotic mass-energy projector
defines the occupied subspace. Conjugating with the y-axis rotation
of angle theta=atan(q/M) gives the exact generator

    omega sigma3 + g0 sigma2,
    g0=-theta_dot/2=q(Mdot+HM)/(2omega^2).

The geometric term is therefore mandatory. The physical energy operator
is the original Hamiltonian, not the extra basis connection.

For fixed p, the difference from the asymptotic Hamiltonian M_asym sigma3
has integrable norm. For T>=max(1,tau), q<=p/|t|^4 and
|M-M_asym|<=Delta tau^8/(8|t|^8), hence the wave-operator error is at most

    p/(3T^3)+Delta tau^8/(56T^7).

The finite-mode wave operators and their adjoints converge in norm.
Their limits are unitary, and bounded convergence over momentum gives
strong one-particle limits on L2. This is not uniform convergence in p
and says nothing about infinite-volume Fock implementability.
Both asymptotic masses remain strictly positive. The original physical
mass rotation and every further fixed finite rotation tend to the
identity at the two time ends.

The geometry and mass source are prescribed here. This does not make
them solutions of the polynomial parent or bound its other sectors.
