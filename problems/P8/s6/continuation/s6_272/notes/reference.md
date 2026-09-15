# Actual near-clock reference and the same full covariance

The S6.270 broad coefficient box was deliberately conservative.
Here use the actual bound |Theta|<=5T, |H|<5T, |E|<=1/2,
|ell|<=1/10, 1<Jc<100, |A|<=10^-390 and
|Tcorr|<3 x 10^-390. Both fixed profiles remain in A and Tcorr.
The physical wave number is q=P^2/a^2; a=(1+u^2)^2 lies in [1,2].

Use the fixed initial scalar balance

    B=diag(sqrt(P),sqrt(P),1/sqrt(P),1/sqrt(P)) T_central(0).

It is symplectic and time independent. Conjugate the complete fixed
canonical generator by B; do not use a moving balance without its
connection. All sixteen actual entry bounds are evaluated. Their sum
is below 10^68. In particular, the terms growing with P are multiplied
by their actual small Theta or profile coefficients, not set to zero.

For tensor shape variables the fixed balance is
diag(sqrt(P)/2,2/sqrt(P)); both shape/dual factors are retained and the
result is the ordinary balanced tensor generator. Each tensor entry
sum is below 10^67. Both tensor polarizations remain.

For the heavy field keep the fixed S6.240 Omega_star, not the physical
instantaneous dispersion. The full block is bounded by

    Omega_hi+(8n+2P^2)/Omega_lo < 10^102.

The physical dispersion still has n+P^2/a^2 and is not replaced by
Omega_star^2=n+P^2/16.

For Proca use the complete fixed initial matrix
B_v=sqrt(I+kk^T/m^2), Omega=sqrt(m^2+P^2), and both zeta normalizations.
In the full longitudinal/transverse basis the balanced upper block is
Omega B_v^-1 K(a) B_v^-1 and is bounded by Omega. Its lower block
B_v V(a) B_v/Omega is bounded by 2Omega on 1<=a<=2. Thus all six entries
sum to at most 9Omega<18P. This keeps the longitudinal sector and its
original Gauss term, not a transverse-only estimate.

Sum six real basis functions, the full scalar block, both tensors,
heavy and all three vector polarizations: the balanced generator norm
is below 10^103.

Let V0=S0 S0^T/2 be the same full pure unit-CCR covariance. Every diagonal
variance in this fixed full balance is bounded by 10^21. This gives

    ||B_full S0||^2 <= 2 trace(V_bal) <= 192 x 10^21.

No cross covariance has been set to zero. Symplecticity gives equality
of the operator norms of B_full S0 and its inverse, so the full whitened
generator is bounded by their product times the balanced generator,
less than 10^128. Tensor and vector orthogonal frame changes do not
alter this norm argument.

For M_w=S0^-1 S_u S0, the full integral equation gives

    ||M_w-I|| <= exp(GT)-1 <= 2GT < 1/50.

Hence the complex initial ball 4R maps strictly into the S0-whitened
image ball 8R on the whole longer interval. This is an exact free-flow
estimate for the unchanged prepared state, not a new vacuum or a
first-order approximation to the interacting propagation.
