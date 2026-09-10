# Positive measure and exact subtracted identity

Use x,z in [0,1], d=xM+(1-x)^2, and

w_L=x^2(1-x),       h_L=(1-x)^2 z(1-z),
w_H=x(1-x)^2,       h_H=x^2 z(1-z).

The already derived loop form factor is
g/(16pi^2) integral sum_i w_i/(d-h_i t).
Its fixed kinetic counterterm subtracts the value at t=0.

On the interior, push forward the positive measure
g w_i dx dz/(16pi^2 h_i) by tau=d/h_i. Call it dnu_i(tau).
Parameter boundaries are treated by limiting integration, not by
dividing by a vanishing h there. The exact resolvent identity gives

F_R(t)-1 = t integral dnu(tau)/[tau(tau-t)],
F_R'(0) = integral dnu(tau)/tau^2,
dnu=dnu_L+dnu_H.

The inverse-first moment equals the finite unsubtracted F_loop(0)
already controlled by S6.112/S6.115. The inverse-second moment is the
finite positive slope. These facts do not require finite total mass
of dnu. On compact sets away from its positive support, the ratio
tau/(tau-t) is uniformly bounded; the inverse-first moment supplies
domination. Tonelli applies to the nonnegative real moments.
These observations justify both the analytic subtracted identity and
the moment interchanges.

The positive gaps

d-4h_L = xM+(1-x)^2(2z-1)^2,
d-4M h_H = Mx(1-x)+(1-x)^2+Mx^2(2z-1)^2

put the supports at [4,infinity) and [4M,infinity). Interior points
approach each endpoint and give positive weight immediately above it.
There is no threshold delta atom. The explicit densities in cuts.md
are locally integrable and vanish continuously at their thresholds.

For T on an open cut, the positive causal prescription gives
Im F_R(T+i0)=pi dnu/dT. The fixed real subtraction has no discontinuity.
This positivity belongs to this particular one-loop vertex measure;
it is not the optical theorem for the complete scattering amplitude.

The two cut integrals reproduce the slope exactly at this retained
order, including arbitrarily high transfer in the mathematical loop
integral. That does not assert a UV-complete gravitational theory or
bound higher loops at all energies.
