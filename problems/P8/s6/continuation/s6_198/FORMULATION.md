# S6.198 formulation

## Actual reference decomposition

Retain the unchanged conditional Gaussian Proca sector,
mass1000, kappa1e800, all-order state, finite prescription
and unit CD clock. Real smooth compact spatial tracefree
D,Gamma may overlap in time; Gamma vanishes near the
common initial Cauchy surface.

Use the actual constant-alpha W8 reference from S186/S197.
For each complete creation pair let Theta'=Omega=W_k+W_l,
g=1/Omega, b=a_Gamma Gammahat and Lb=(g b)'.
Six exact retarded integrations give all6 boundary terms
and a sixth-order bulk, with no initial endpoint left over.

The j5 endpoint and sixth bulk have a joint full-continuum
limit F. Define S6[Gamma]^2=sum_(j=0)^6||partial_t^j Gamma||L2^2.
For nonzero tests,

    |F|<1e48 ||D||L2 S6[Gamma];
    |F-F_K|<1e52 ||D||L2 S6[Gamma]/K, K>=1000.

There is no spatial Fourier cutoff. The sixth inverse-frequency
power is integrable; the analogous fifth-step absolute bulk
majorant is not. No cancellation of the exact fifth bulk is
ruled out by that majorant observation.

## Complete known finite piece

At every common regulator the actual current is exactly

    J_actual,K=R_state,K+C_ref,K+sum_(j=0)^4 B_j,K+F_K.

R_state is the complete S197 actual-minus-reference memory
and contact. C_ref is the FULL reference metric contact,
and every B_j is the explicit equal-time term.

With M[D]^2=||D||L2^2+||grad D||L2^2 and
N61[Gamma]^2=S6[Gamma]^2+||grad Gamma||L2^2, the known
finite piece R_state+F has bound2e48 M[D]N61[Gamma]
and regulator error2e52 M[D]N61[Gamma]/K.
Both canonical tensor factors give8e-752 and8e-748/K.
All quantities vanish exactly if a test is zero.

## Boundary

The full contact plus five earlier endpoints still need the
ORIGINAL covariant spatial matching. Equal-time source jets
do not establish spatial locality, a reference Ward identity
or the cancellation of any odd endpoint. The contact and
memory retain their different regulator regions.
No counterterm is reset. No full mixed inverse, nonlinear
background, stability, physical cutoff or V/G/B closure is
established. Original P8 remains OPEN.
