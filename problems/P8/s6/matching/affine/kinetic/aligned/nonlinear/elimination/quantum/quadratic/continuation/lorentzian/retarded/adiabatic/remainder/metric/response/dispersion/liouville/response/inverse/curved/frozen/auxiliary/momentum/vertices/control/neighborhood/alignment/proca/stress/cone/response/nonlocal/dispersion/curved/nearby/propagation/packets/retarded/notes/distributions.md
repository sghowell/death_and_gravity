# Exact retarded kernel: two shells plus an L2 remainder

Use the convolution-kernel Fourier convention

    G(x)=(2pi)^(-3) integral exp(i k_vector.x) G_hat(k_vector) d^3k.

For any S>0,

    inverse_Fourier[sin(kS)/k] = delta(|x|-S)/(4pi|x|).

One can verify the normalization directly by the forward
transform. Angular integration is 2pi integral_-1^1
exp(-ikS mu)dmu. The radial delta integral leaves S/(4pi).
Their product is sin(kS)/k. remainder.py checks this angular
identity and its finite k=0 limit. The standard retarded
wave-kernel Fourier construction is also given in
[Manchester's wave-equation notes](https://oer.physics.manchester.ac.uk/AM/Notes/jsmath/Notesse12.html).
The variable-time shell radii and amplitudes here are derived
from the actual system, not imported from a constant-coefficient
wave equation.

## High-frequency remainder

Put C=10^44 and K0=4*10^31. S6.91's actual two-time error gives

    integral_{k>=K0} d^3k |G_hat-G_hat^0|^2
       <=4pi C^2 integral_K0^infinity k^-2 dk
       =4pi C^2/K0.

Plancherel therefore bounds the squared spatial L2 norm of this
high-frequency remainder by C^2/(2pi^2 K0).
This only uses the absolute multiplier bound. No unproved
frequency-derivative symbol estimates or smooth-error
assumption are needed.

## The entire low-frequency ball is included

For 0<=k<=K0, q<=4K0^2 since R>1/2. Apply the actual parent's
coefficient bounds to every coefficient of the polynomial
original-chart matrix M_X. remainder.py records the complete
nonnegative entry majorants and their exact maximum row sum B.
It is finite, though very large. Gronwall gives

    ||T_X(t,s;k)||<=exp(T B).

The source factor is below 16; the original output row has
norm one. Since |sin(kS_j)/k|<=S_j<4T and each A_j<10^6,

    |G_hat-G_hat^0| <= M_low
       :=16 exp(T B)+8*10^6 T.

The compact-ball squared spatial L2 bound is consequently
K0^3 M_low^2/(6pi^2). It is a finite existence bound, not a
numerically useful low-frequency stability estimate.
No frequency range, finite-q Hessian zero or contact term
has been removed from the exact classical response.

Combining the disjoint frequency pieces proves

    G(t,s;x)=A_c delta(r-S_c)/(4pi r)
             +A_m delta(r-S_m)/(4pi r)+R(t,s;x),
    R(t,s;.) in L2(R^3).

The L2 bound is uniform on ordered pairs in I. The original
ODE and the explicit leading transfer are continuous in both
times at every frequency. The same integrable majorants give
L2 continuity of R by dominated convergence. In particular,
R is locally L2 jointly in t,s and x, away from no additional
spatial exception. Translating x to x-y preserves this property
on bounded product source/detector neighborhoods.

## Why the front cannot be cancelled

For s<t, both radii are positive and S_c>S_m. Near a point of
the clock shell, the matter shell is absent. Its defining
normal coordinate z=r-S_c is smooth with nonzero radial
derivative, and A_c/(4pi r) is smooth and positive. The local
kernel has the form a times delta(z), plus a locally L2
function, with a nonzero on an open tangential patch.

A nonzero surface delta is not locally L2. To see this without
a microlocal theorem, choose a fixed nonnegative tangential
bump and a smooth normal bump psi(z/epsilon) with psi(0)=1.
Its pairing with the surface delta is independent of epsilon
and nonzero. Its ordinary L2 norm is O(sqrt(epsilon)).
Cauchy-Schwarz makes its pairing with any locally L2 function
tend to zero. Thus the L2 remainder cannot cancel the clock
shell on any neighborhood of that patch.

This is a written distribution proof with exact normalization,
coefficient and integrability anchors, not proof-assistant
formalization. It establishes support outside the matter cone,
not that the entire remainder is smooth or that every possible
lower-order singularity has been classified.
