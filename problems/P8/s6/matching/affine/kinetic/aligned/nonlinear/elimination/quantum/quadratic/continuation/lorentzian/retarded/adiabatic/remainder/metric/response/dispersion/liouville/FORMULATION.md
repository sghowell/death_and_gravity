# P8 S6.70: exact acoustic-time vector reduction

This checkpoint leaves the original clock, physical matter metric,
selected all-order state and fixed subtraction prescription unchanged.
Acoustic time is an auxiliary mode coordinate, not a different
physical frame or a new causal-cone claim.

## Statement

For smooth positive physical N,a,am,bm and positive comoving k,
the two transverse and the constrained longitudinal oscillators
admit exact positive Liouville changes of time and canonical
normalization. They preserve the Wronskian and the physical readout.

The transverse equation has potential a^2 m^2 bm in its conformal
time. The longitudinal acoustic time has rate N sqrt(bm/am)/a.
Its pump is proportional to

    A k/sqrt(k^2+U), A=a(am bm)^(1/4), U=a^2 m^2 am.

The exact longitudinal potential is

    V_k=U-A''/A+(b U'+U''/2)/(k^2+U)
                      -3 U'^2/[4(k^2+U)^2], b=A'/A,

where primes here mean longitudinal acoustic derivatives.

On the original clock and I=[-1/2,1/2], the momentum remainder R_k
is nonnegative and bounded by 88 m^2/(k^2+m^2). Its scalar principal
potential is at least m^2-275/16, hence positive at the selected
m=1000. The potential, though not the k-dependent canonical map,
has a regular k=0 limit reproducing the transverse potential.

For smooth physical lapse/log-scale sources zero on an initial
neighborhood, the potential-remainder first variation at fixed
acoustic time includes the original retarded coordinate shift.
The joint physical source C2-to-potential-C0 bound is

    ||delta R_k|| <=18000 m^2/(k^2+m^2) ||(n,zeta)||_C2.

The source chart chi=zeta+(alpha+beta)n/4 and
r_sigma=[1+(beta-alpha)/2]n-zeta is regular through the bounce.
The second derivative of chi gives exactly the rank-one source
direction already measured by S6.69's ultraviolet response pole.

## Boundary

This is a potential-level reduction, not a renormalized stress
estimate or an exact replacement by one minimally coupled scalar.
The momentum-dependent pump and both distinct acoustic times remain.
The subtraction/readout variation, complete causal matrix inverse,
spatial response, independent initial states, quantum stability/cones,
interactions, other loops, cutoff, finite Wilson matching and V/G/B
are not supplied here. Original P8 remains open.

Evidence is exact symbolic/rational computation and written continuous
mode/potential arguments, not proof-assistant formalization.
