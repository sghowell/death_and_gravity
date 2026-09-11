# Exact frames, boundary states and the complete transition

The construction uses exact two-level changes of basis.
The general iterative method appears in
[Barbero et al., section IV, equations4.6-4.9](https://arxiv.org/pdf/1805.05107).
Our constants and remainder bound are derived here.

After the physical mass rotation, in x=t/tau units,

    H0=e0 sigma3+g0 sigma2,
    e0=tau sqrt(p^2+M^2),
    g0=pDelta s'/(2(p^2+M^2)).

For the opposite mass profile or helicity the corresponding
sign changes in g0; the absolute estimates are identical.

For H=e sigma3+g sigma1, the y-axis rotation that
diagonalizes H induces connection
-minus one half the derivative of atan(g/e), on sigma2.
For H=e sigma3+g sigma2, the corresponding x-axis
rotation induces the positive half derivative, on sigma1.
Both follow by computing U^dagger H U-iU^dagger U'.

Consequently the exact four-step recursion is

    e_(j+1)=sqrt(e_j^2+g_j^2),
    g_(j+1)=(-1)^j partial_x atan(g_j/e_j)/2,
      j=0,1,2,3.

The off-diagonal axes alternate 2,1,2,1,2.
No connection term is omitted. The derivative may also
be evaluated as (e_j g_j'-g_j e_j')/[2(e_j^2+g_j^2)].

The complex-disk proof shows that all g_j vanish at both
infinities, while the real frequencies remain positive.
Each additional rotation therefore tends to the identity.
The original physical in/out eigenbases and transition
probability are preserved. The initial mass rotation
itself is not falsely set to the identity.

Remove the final diagonal phase, leaving the exact
off-diagonal generator with norm |g4(x)|. Unitarity
and its integral equation imply that its full
off-diagonal transition is bounded by integral |g4|.
This retains every transition-coupling Dyson order;
it is not an assertion about all Feynman loops.

Using integral w=2 gives

    |beta_p|<=integral |g4|dx
       <=4*1024^4 pDelta/(tau^4 E^6)
        =2^42 pDelta/(tau^4 E^6).

This falls as p^-5 and vanishes at p=0. Mode unitarity
also bounds it by one. The convenient uniform allowance

    |beta_p|<=2^42 Delta/(2tau^4 m0^5)

uses p/E^2<=1/(2m0) and E>=m0. The actual value is
approximately 6.93705445956e-391, below 1e-390.

Independent tests check both rotations and their signs,
evaluate the exact recursive derivatives using finite
Taylor jets, compare with direct differentiation, and
test the complex disk bounds. Finite jets here evaluate
derivatives of the full functions; they do not truncate
the physical Hamiltonian or evolution.
