# Dimensional insertion before the outer MS limit

Use d=4-2epsilon with trace(1_Dirac)=4 and
the common MS-bar loop normalization
e^(gamma epsilon) mu^(2epsilon)/Q.
At the end mu=mF. A scalar bubble is

 B_D(m,m;s)=e^(gamma epsilon) Gamma(epsilon)
   mu^(2epsilon)/Q
   integral_0^1 [m^2-x(1-x)s]^(-epsilon) dx.

The literal Dirac trace and a numerator
dot-product identity give

 f_D(s)=2NY[(4m^2-s)B_D(m,m;s)-2 Tad_D(m)].

At mu=m the zero-momentum bubble has pole
1/epsilon and finite part zero, while the
tadpole relation Gamma(epsilon-1)
=Gamma(epsilon)/(epsilon-1) gives the frozen
mass-inverse pole 6m^2 and finite part 2m^2
after stripping 2NY/Q. Thus this normalization
agrees with S6.131 and the shared MS reference.

Keep the exact D-dimensional mass/residue
subtraction, not only its four-dimensional
finite value. Its tadpole cancels:

 f_R,D(s)/(2NY)
 =(4m^2-s)[B_D(s)-B_D(1)]
  -(s-1)(4m^2-1) B_D'(1).

This is the same physical mass-one reference
used in the frozen insertion. It removes its
whole affine proper term before any outer
limit, retaining its dimensional continuation.

For T=4m^2, beta=sqrt(1-T/v), set
x=(1+beta u)/2 on the cut. Then
v x(1-x)-m^2=v beta^2(1-u^2)/4.
The discontinuity and the Euler beta integral,
with Gamma(epsilon) Gamma(1-epsilon)
sin(pi epsilon)/pi=1, give

 rho_D(v)/C
 =e^(gamma epsilon) mu^(2epsilon)4^epsilon
  sqrt(pi)/(2 Gamma(3/2-epsilon))
  v^(1-epsilon)(1-T/v)^(3/2-epsilon),

where C=2NY/Q. At epsilon=0 this is exactly
the frozen positive scalar insertion density.

For real 0<epsilon<1 the two-point analytic
function with its affine terms removed has
the convergent Stieltjes representation

 -f_R,D(s)/(s-1)^2
   =integral_T^infinity rho_D(v)/
                      [(v-1)^2(v-s)] dv.

The bubble parameter formula supplies the
cut and large-complex-energy decay of this
one-loop function; no unproved high-energy
property of the complete theory is invoked.
This is a regulator identity, not a new
physical spectrum or normalized exact measure.

The outer zero-momentum reference is therefore

 F_D=integral_T^infinity
       rho_D(v)/(v-1)^2 B_D(v,1;0) dv.

It includes both loop factors and converges
in the stated positive regulator range.
All finite epsilon-times-pole products are
retained before taking its MS finite part.

The mu powers also match the parent vertices:
each C(z)=-L+g/(M-z) has bare mu^(2epsilon),
as does Y; after factoring the external
quartic mu^(2epsilon), two loop factors
mu^(2epsilon) remain. The same counting
gives the local cubic and heavy-mass
counterterms in the frozen parent map.
