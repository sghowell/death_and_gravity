# Local nonlinear constraints of the new retained action

Use the timelike scalar-unitary chart, with the exact time and
spatial boundaries in notes/chart.md. The physical lapse N=1/s,
normal vector component T_n, shift, and their momenta must all
be retained. Maxwell's Legendre transformation and spatial
integration by parts give the normal-vector Gauss multiplier
N*T_n*div(pi); no lapse derivative is introduced. The free M1
field remains a field during metric variation.

Let v denote hat-volume ratio R^(-3/4), delta=R-1 and
c=-3H s delta. Apart from the shear, vector electric and matter
blocks, the velocity density has the exact form

    L=a K^2+b K+f+v(T_n-delta K-c)^2/2,
    a=-v R/3.

Here b and f include BOTH contributions from the primitive I.
The effective trace coefficient a+v delta^2/2 is nonzero because

    gamma_eff=1-3(R-1)^2/(2R),
    gamma_eff-1/4=3(R-1/2)(2-R)/(2R)>0.

Thus on the whole coefficient domain the ten-velocity block is
nonsingular: one unreduced negative gravitational trace direction,
five positive shears, three positive electric directions and the
positive matter direction. Unreduced trace inertia is not a
physical ghost verdict.

With trace momentum p and Gauss density j (volume normalization
as in the code), exact joint elimination gives

    K=(p-b-delta j)/(2a),
    T_n=delta K+c-j/v,
    H=(p-b-delta j)^2/(4a)-f+j^2/(2v)-c j.

This keeps the Gauss term; dropping it changes the system.
The temporal secondary Jacobian is nonzero by gamma_eff>1/4.
The lapse equation remains a separate condition.

On the actual clock, the new S and its first variation vanish.
The joint lapse/normal-vector Hessian at the background therefore
equals diag(-2J,-1), with all physical kappa/volume factors positive.
J is computed afresh from the source-pinned target F clock jets,
original chi momentum, and full ADM boundary:

    J=P(u^2)/[800(1+u^2)^18].

The polynomial P has degree17 and strictly positive coefficients,
recorded exactly in the report; its constant is1215 and leading
coefficient3200. Hence J>0 at EVERY finite u, including u=0.
No division by H or gamma-crossing Theta has occurred.
The required primitive jets are I_N=0 and I_NN=-3h_u/h^3;
omitting the latter changes the answer.

Continuity and the implicit function theorem give a nonempty
local canonical finite-jet neighborhood around each finite
clock point where the two auxiliary equations solve for lapse
and normal vector, with their preservation determining the
corresponding multipliers. The exact spatial diffeomorphism
constraints persist. In the unitary chart there are15
configuration variables:6 metric,4 vector,1 matter,1 lapse,
3 shifts. Six first-class constraints (shift momenta and spatial
diffeomorphisms) and four second-class constraints (the two
auxiliary primary/secondary pairs) leave14 phase-space dimensions,
or seven physical modes.

This is a local nonlinear constraint count. J tends to zero in
the tails, so no uniform neighborhood radius follows from its
strict finite-time positivity. No global nonlinear solution,
stability, controlled cone neighborhood or UV cutoff is inferred.
The vacuum is checked quadratically in covariant variables,
not by extending scalar-unitary gauge through du=0.
