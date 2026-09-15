# Exact massless endpoint and nonuniform forward limit

For beta=i*b, b>0, the COMPLETE forward angular integrand is

    [(1+b^2)^4+b^8(1-x^2)^4]/(1+b^2 x^2)^2.

Its exact polynomial division and two elementary primitives on x in[0,1]
give the full shape

    F_b=4+31b^2/3+118b^4/15+b^6
        +(b^7-6b^3-8b-3/b)atan(b).

Use w=sqrt(s/(4mu)) and atan(sqrt(1-w^2)/w)=pi/2-asin(w).
The analytic power series near w=0 proves

    kappa^2 rho/mu^2 =
      1/(32w^3)-7/(64w)+41/(80pi)
      -13w/256-w^2/(60pi)+O(w^3).

In particular the unweighted low-cut density is not integrable at0.
The remainder assertion follows from the convergent analytic asin and
square-root expansions on every sufficiently small complex disk with
their fixed branches, not from numerical sampling.

For strict |z|<1 put h=-tau^2. The two negative-real seeds have constant
terms -acos(z)/sqrt(1-z^2) and -(pi-acos(z))/sqrt(1-z^2).
Their exact first derivatives in tau^2 are, respectively,

    1/(1+z)+acos(z)/[(1+z)sqrt(1-z^2)],
    1/(1-z)+(pi-acos(z))/[(1-z)sqrt(1-z^2)].

The denominator and atan argument are analytic in tau^2 uniformly on
every compact strict-angle set. Together with L0=-pi/(2tau)+1+O(tau^2),
the complete angular formula gives the coefficients in FORMULATION.md.
Here s=4mu*tau^2/(1+tau^2). The O(1) remainder is uniform on those compact
strict-angle sets, not at z=+/-1. At z^2=1/2 the1/s coefficient vanishes;
this does not contradict the formula or imply a generic cancellation.

At fixed0<t<4mu, z tends to1-t/(2mu). The complementary-kernel identity

    M(t)+M(4mu-t)=pi/[2sqrt(t(4mu-t))]

also derives the leading1/s coefficient directly in invariant variables.
At t=0 the endpoints merge and the stronger s^(-3/2) forward term results.
Consequently neither the forward limit nor an uncompensated endpoint
integral can be interchanged with fixed-angle asymptotics.
