# General Ward identity and exact soft normalization

An independent invariant-index algebra derives the gauge contraction without
choosing a physical kinematic point. Use basis p1,p2,p3,k,xi with
p1^2=p2^2=p3^2=mu,k^2=0,p1.p2=a,p1.p3=b,pi.k=ui.
Set p4=-p1-p2-p3-k and
p2.p3=-mu-a-b-u1-u2-u3, which enforces p4^2=mu.
The pi.xi and k.xi are arbitrary. No four-dimensional Gram determinant
constraint is needed; an identity in this larger invariant algebra in
particular holds on the physical four-dimensional locus.

Represent each tensor term by metric edges and vector endpoints. Every
contracted index occurs twice. Its connected components are paths joining
two vectors (giving their Minkowski dot product), or closed metric loops
(giving dimension4). This independent contraction algorithm expands the
same canonical vertices as the component implementation.

Replace the external polarization by k_lower xi_lower+xi_lower k_lower.
The emitted scalar vertex is exactly -4(pi.k)(pi.xi); division by its scalar
propagator2pi.k gives-2pi.xi. For the ij|lm partition multiply the sum of
the four external graphs, two seagulls and cubic graph by K_ij^2 K_lm^2.
The entire resulting polynomial in mu,a,b,u1,u2,u3 and the four gauge
contractions vanishes. No point samples are used in this cancellation.
Bose relabeling proves the other two partitions. Reversing the cubic sign
gives a nonzero polynomial, providing an explicit negative control.

A separate invariant expression for epsilon=xi_lower xi_lower is compared
with all seven component diagrams at each of three rational configurations
and each channel. This checks the full vertex, not merely its gauge part.

For physical TT epsilon, k.epsilon=0 and tr(epsilon)=0, so
V3(epsilon;pi,-pi-k)/(2pi.k)=-pi.epsilon.pi/(pi.k)=-J_i.
The remaining graviton-exchange vertex product is minus the corresponding
four-point amplitude. Hence, exactly,
sqrt(kappa)^3 M_G5=sum_i J_i G_i+seagulls+cubic,
where G_i is the whole off-shell-shifted four-point gravity numerator sum.
At k->0 and fixed nonzero hard transfer, G_i->G0=kappa A_G.
Thus M_G5->A_G sum_i J_i/sqrt(kappa).
The next note bounds this difference rather than inferring uniformity
from the limiting identity.
