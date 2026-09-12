# Uniform coefficient bound from negative-axis Cauchy discs

Fix real p, q=|p|^2, and the closed complex disc |z+q|<=m^2.
For s>=4m^2,

    |s-z|>=3s/4,
    |z|<=q+m^2,
    ||A(z)||op<=|z|+q<=2(q+m^2).

For a complex symmetric matrix D, Q0 acts by A tr(A D)/3 and Q2 by
A D A^T-A tr(A D)/3. Frobenius Cauchy-Schwarz and
||A||F<=sqrt(3)||A||op give

    ||Q0||<=||A||op^2,
    ||Q2||<=2||A||op^2.

No Euclidean orthogonal-projector norm is assumed at nonzero p.
Consequently the analytic matrix-valued F satisfies

    ||F|| <= (32/3)(I2+I0)(q+m^2)^3
           =(q+m^2)^3/(80pi^2 m^2),

because I2+I0=3/(2560pi^2m^2). These integrals and bounds
are uniform on the disc, so the operator-valued Cauchy formula
applies. Its Taylor coefficients about z=-q obey

    ||A_r(p)|| <=(q+m^2)^3/(80pi^2 m^(2r+2)), r=0,1,2.

The same argument applied only to s>=Lambda uses
I2_tail+I0_tail<=1/(96pi^2 Lambda), giving

    ||A_r-A_r,Lambda|| <=(q+m^2)^3/(9pi^2 m^(2r) Lambda).

This is a controlled improper tail, not a finite sampling estimate.
It holds at every real p.

In the weighted norm

    Y^2=sum_(r=0)^2 m^(-4r)||(m^2-Delta)^3 partial_t^(2r)Gamma||L2^2,

Cauchy-Schwarz across the three time orders yields the conversion
bound sqrt(3)/(80pi^2m^2)||D||L2 Y and tail
sqrt(3)/(9pi^2 Lambda)||D||L2 Y.

For m>=1, (q+m^2)^3<=m^6(1+q)^3 and m^(-2r)<=1.
Thus at m=1000 the standard Z norm in FORMULATION.md gives the
conservative displays 1e10 and 1e17/Lambda. This norm retains six
spatial and four time source derivatives; its derivative loss is
not hidden in the numerical constants.
