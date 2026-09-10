# Exact constant-background one-loop potential remainder

Keep the same canonical polynomial model and fixed parameters.
For constant Phi, the classical heavy solution is
H_*=-G Phi^2/(2mu^2). In Euclidean momentum squared y>=0,
the full two-field Hessian, differentiated before this substitution,
has heavy block y+mu^2 and reduced light Schur block

    y+1+F(y) Phi^2,
    F(y)=A-G^2/(y+mu^2),
    A=(lambda4-G^2/mu^2)/2,
    delta=lambda4-3G^2/mu^2.

It follows exactly that

    F(y)=delta/2+G^2 y/[mu^2(y+mu^2)],
    delta/2<=F(y)<A, F'(y)>0.

Thus the full constant-background Hessian is positive for every
Phi and y. This does not state positivity for all inhomogeneous
backgrounds. The heavy determinant is field-independent.
The Schur kernel still contains the full heavy propagator.

At one loop the stationary light effective potential is obtained
by the full Hessian on the classical heavy solution. A one-loop
shift of that solution affects its value only at higher loop order,
since the tree action is stationary there. This also follows from
the exact finite-regulator light functional.

Subtract the vacuum constant and the Phi^2 and Phi^4 Taylor
coefficients at Phi=0, once for this new flat model. With

    z=F(y) Phi^2/(y+1), R(z)=log(1+z)-z+z^2/2,

the resulting remainder is V1_ren=(1/2) integral d^4q/(2pi)^4 R(z).
These are fixed effective-potential renormalization conditions,
not on-shell pole/residue conditions or a selection of b2.
They are unrelated to frozen P8 state or tadpole prescriptions.

For z>=0,

    R(z)=integral_0^z t^2/(1+t) dt,
    z^3/[3(1+z)]<=R(z)<=z^3/3.

Use F<=A, z<=A Phi^2 and F>=delta/2. The elementary
four-dimensional radial integral is

    integral d^4q/(2pi)^4 (q^2+1)^-3
       =(1/(16pi^2)) integral_0^infinity y/(y+1)^3 dy
       =1/(32pi^2).

It proves convergence after the stated local subtractions and gives

    delta^3 Phi^6/[1536pi^2(1+A Phi^2)]
       <=V1_ren<=A^3 Phi^6/(192pi^2).

No expansion of a propagator under an unbounded integral was used.
The bounds hold for every finite real Phi; all equal zero at Phi=0,
and the lower bound is strictly positive elsewhere.

For an elementary rational pi lower bound, integrate the finite
identity 1/(1+t^2)=sum_(j=0)^7(-t^2)^j+t^16/(1+t^2)
on 0<=t<=1. Four times the polynomial integral is strictly
greater than 3, so pi>3 without a floating-point estimate.
At the exact selected couplings, A<10^-205 and the sharper
actual A^3/1728<10^-618. The last ceiling is a rational
calibration, not a rounded decimal assertion.

A stronger ceiling preserves the heavy-mass dependence:
write F(y)=a+b y/(y+mu^2), a=delta/2, b=G^2/mu^2.
The identity 4(a^3+b^3)-(a+b)^3=3(a-b)^2(a+b)>=0 and
y^4/(y+1)^3<=y give

    integral_0^infinity y F(y)^3/(y+1)^3 dy
       <=4[a^3/2+b^3/(2mu^2)].

The second term uses the exact scaled integral
integral y/(y+mu^2)^3 dy=1/(2mu^2). Consequently

    V1_ren <= [a^3+b^3/mu^2] Phi^6/(48pi^2)
            <[a^3+b^3/mu^2] Phi^6/432
            <10^-815 Phi^6 for Phi != 0.

All inequalities used positive integrands and kept the full
heavy propagator. The improvement is an actual mass-weighted
loop estimate, not a formal exchange of a large-mass series
and an unbounded momentum integral.

This positive, bounded one-loop remainder cannot overturn the
constant-field tree global minimum in this fixed subtraction scheme.
That statement is about the tree-plus-one-loop potential only.
There is no uniform relative higher-loop estimate for arbitrary
large Phi, no nonperturbative vacuum theorem, no quantum pole-mass
certificate and no on-shell four-point/derivative bound.
