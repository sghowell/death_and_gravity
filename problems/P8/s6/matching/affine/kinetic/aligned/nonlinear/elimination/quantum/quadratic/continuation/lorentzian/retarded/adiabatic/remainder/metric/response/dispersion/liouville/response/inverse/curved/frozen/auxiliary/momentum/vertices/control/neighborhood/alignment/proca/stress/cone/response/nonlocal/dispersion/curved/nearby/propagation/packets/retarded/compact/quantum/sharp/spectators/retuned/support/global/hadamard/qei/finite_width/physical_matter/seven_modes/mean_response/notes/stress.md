# Conserved physical stress and the spatial-only lapse chain

Let Z=(W_i,Pi_i), and G=diag(B_P,A_P) be the actual Cartesian
Proca Hamiltonian Hessian of S6.103, so H_P=Z^T G Z/2.
At fixed canonical fields define

    R=G/a^3,  P=-G_a/(3a^2),
    r=(1/2) integral tr(R Delta C) d^3k/(2pi)^3,
    s=(1/2) integral tr(P Delta C) d^3k/(2pi)^3.

These are the physical normal density and one-third spatial
stress trace. The pressure formula follows from variation of
the physical spatial scale, not from a guessed equation of state.

There are four positive density pieces: electric momentum,
temporal/longitudinal momentum, spatial vector mass, and
magnetic curl. Their pressure weights are respectively
1/3, 1, -1/3 and 1/3. Positivity of Delta C therefore gives
r>=0 and -r/3<=s<=r at every time. The magnetic matrix is
the square of the Cartesian curl map; no angular frame
or four-scalar Proca replacement is used.

The native matrix identity

    a H R_a+M_P^T R+R M_P+3H(R+P)=0

proves r'+3H(r+s)=0 after differentiation under the finite
momentum integral. The source is thus actually conserved,
not assigned a pressure to force formal conservation.

The inherited transformation is SPATIAL:
h_physical=e^2 h_hat with the same N and shift.
Holding W,Pi and hat_a fixed, the physical Hamiltonian is

    N H_P(e(N,u) hat_a).

Since e_N|N=1=1/(2h), its N derivative per background
hat volume is r-3s/(2h). A false four-dimensional
conformal lapse replacement N->N e would instead add
r/(2h). The code retains that nonzero difference as a
negative control and checks the correct full source.

At the center the induced scalar matter density is
-2(r-3s/2)/405. All 36 entries of its six-channel Hessian
match S6.102's independent constrained calculation.
Adding the direct vector density gives

    delta rho_total=(403r+3s)/405 >=134r/135.

This comparison holds for the declared zero-anchor mean
scale/trace data. It does not make the isolated scalar
second variation positive, nor establish a general QEI.
