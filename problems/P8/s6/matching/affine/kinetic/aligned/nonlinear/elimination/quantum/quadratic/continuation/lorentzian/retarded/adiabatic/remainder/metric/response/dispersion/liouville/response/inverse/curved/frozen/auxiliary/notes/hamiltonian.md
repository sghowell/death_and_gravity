# Full local Hamiltonian and its holomorphic coefficient bound

Use the positive spatial hat metric of S6.45. Work in the
existing dimensionless clock units and divide the normal
Hamiltonian density by sqrt(det hat_h). Canonical momenta
are held fixed in lapse derivatives.

The trace variable is p=2*pi_trace/(3*sqrt(det hat_h)),
with background p0=-2H. The matter momentum per hat volume
has background ell. Define nine independent invariant
coordinates, used only to describe this Hamiltonian:

1. dp=p+2H;
2. dc=p_chi-ell;
3. j=D_i Pi^i/sqrt(det hat_h);
4. S=pi_TF,ij*pi_TF^ij/det(hat_h);
5. e=hat_h_ij*Pi^i*Pi^j/[zeta*det(hat_h)];
6. bmag=zeta*F_ij*F_kl*hat_h^ik*hat_h^jl;
7. v=hat_h^ij*W_i*W_j;
8. cgrad=hat_h^ij*D_i chi*D_j chi;
9. R=R3_hat.

Their common background value is zero. For real physical
data the squares have their usual positivity, but allowing
independent complex values gives a larger polydisc on which
Cauchy estimates are available. The metric remains positive
when interpreting real data. No spatial momentum constraint
is discarded or solved by this coordinate description.

The e and bmag definitions retain inverse-zeta electric
normalization and zeta magnetic normalization. They are
not bounds on unnormalized coordinate components of Pi
and W at arbitrary momentum.

## Actual coefficients

Put s=1/N, h=(1+u^2)^3, r=(h-1+s^2)/h,
eomega=r^(-1/4), U=r^(-3/4), B=-r/2, a=2UB/3.
Let I(u,s) be the original boundary primitive with
I_s=U*12*s^2*f_phi*omega_X and I(u,1)=0.
Then the exact S6.45 coefficients, including the scalar
margin, are

b=U*Blinear-I,
f0=U*(Fnew+epsilon*(1-s^2)^2/h^2)-s*I_phi,
delta=(s^2-1)/h,
c=-3H*s*delta+(3/2)*s*Q.

Q is the original ODE solution, not a free function.
The retained mass functions, for p_affine=sqrt(r)/2, are

gamma_t=9*(2p_affine^3-1)/(22p_affine^3+3p_affine-11),
gamma_s=9*(8p_affine+5)/(-72p_affine^2+88p_affine+55).

With p=-2H+dp and l=ell+dc, the entire post-temporal
Hamiltonian per hat volume is

Hcal=N*[
 (p-b-delta*j)^2/(4a)-f0+gamma_t*j^2/(2U)-c*j
 -S/(UB)+e/(2eomega)+bmag/(4eomega)
 +eomega*v/(2gamma_s)+l^2/(2U)
 +eomega*cgrad/2-eomega*B4*R ], B4=-B.

The exact trace and temporal reconstructions are

K=(p-b-delta*j)/(2a),
T=delta*K+c-gamma_t*j/U.

The code checks these against the original joint Legendre
transformation and checks both Maxwell normalizations.
The shear normalization uses pi_TF, not twice pi_TF.
The curvature and matter terms are the complete S6.45
spatial terms. All these variables are local canonical
jets; none introduces a derivative acting on N or T.

The original clock equation is Hcal_N(1,0)=0. Its lapse
derivative is d0=Hcal_NN(1,0)=-2J_e,
J_e=J+4*epsilon/h^2. The literal boundary primitive is
included in this statement. S6.45's exact polynomial
bound gives J>1/40 on I, so |d0^-1|<20. The margin only
strengthens this inequality.

## Complex domain and the entire scalar function

For every real u0 in I, allow |u-u0|<=1/50 and
|N-1|<=1/100. Then

|u|<=13/25, |1+u^2|>99/100,
|s-1|<=1/99,
|s^2-1|<=199/9801<1/40.

The first lower bound follows from
Re(1+u^2)>=1-(Im u)^2>=1-1/2500.
Thus |h|>(99/100)^3, and for |x+1|<=1/20,
|h-1-x|>9/10. Also |h'|<6.
The ratio r lies within 3/100 of 1. Principal powers
therefore have one consistent holomorphic branch,
and |eomega|,|eomega^-1|,|U|,|U^-1|<2.
Since UB=-r^(1/4)/2, |UB|>1/4 and |(4a)^-1|<2.

The principal square root obeys
|p_affine-1/2|<=|r-1|/2<1/50:
Re sqrt(r)>0 implies |sqrt(r)+1|>=1.
Expanding each mass denominator around p_affine=1/2
gives the exact rational bounds in the report:
gamma_t, gamma_s and 1/gamma_s all have magnitude below 2.
The temporal denominator is bounded away from zero.

The entire frozen scalar F has 33 monomials in its
numerator and denominator 200*(1+u^2)^12. Its coefficient
absolute-value sum at |u|=13/25, |x|=21/20, divided by
200*(99/100)^12, is

101704307704782433883533238264601170673664 /
422661243303360596017265796661376953125 < 1000.

This uses the full covariant scalar, not just its clock
Taylor polynomial. The actual transformed coefficients
also obey |omega_phi|<1/10 and |f_phi|<1/10, since
their factors are respectively
h'*(s^2-1)/(4*h*(h-1+s^2)) and
h'*(s^2-1)/(2*h^2).
Consequently Fnew plus the small margin is below 1002.

## Original source and boundary primitive, not formal jets

For y=x+1 the actual source equation is

Q_y+A Q=Fq, Q(u,0)=0,
A=1/[2(y-1)]+3/[4(h-y)],
Fq=-3h'*y/[4h^2*(h-y)].

On the complex y disc of radius 1/20, |A|<3/2 and
|Fq|<7|y|. Integrating on the straight segment and using
the holomorphic variation-of-constants solution gives

|Q| <= (7/2)*exp(3/40)*|y|^2
     <= (140/37)*|y|^2 < 4|y|^2.

The exponential bound follows termwise from exp(t)<=1/(1-t)
for 0<=t<1. The coefficient functions have no pole in the
joint complex domain, so this is also holomorphic in u.

The displayed exact I_s has magnitude below 2 along the
straight segment from 1 to s. Thus |I|<=2/99<1/40.
For fixed complex N, Cauchy's formula on the u circle
of radius 1/50 gives |I_phi|<=100/99<2 at real u0.
This is not differentiation of an absolute-value bound.
It also proves that I_phi is holomorphic in N.

At real u in I, these bounds give |b|<1, |c|<2,
|f0|<3000 and |K|<27 on max|Y_i|<=1.
For b, the exact Blinear=3s*u*(s^2-1)/(1+u^2)^4
and the I bound suffice. For c use |H|<3 and |Q|<4|y|^2.
For f0 use U<2, Fnew+margin<1002 and |s I_phi|<4.

Finally, the eleven Hamiltonian term bounds before its
outer N are respectively

162, 3000, 2, 2, 4, 1, 1/2, 2, 4, 1, 2.

Since |N|<2, their sum gives |Hcal|<6361<10000,
uniformly for real u in I, |N-1|<=1/100 and
max|Y_i|<=1. This bound includes every spatial term,
the actual source and boundary primitive.

## Independent closed bounce anchor

At u=0 the original Q and I vanish for every nearby N,
but I_phi does not. Direct differentiation of the actual
primitive integrand and its fixed basepoint gives

I_phi(0,s)=-6s^(3/2)-18s^(-1/2)+24.

The code checks its s derivative against the original
I_s,u and checks its value at s=1. Substitution into the
full Hamiltonian gives an elementary expression A(N)+L(N)+Q2(N),
retaining every invariant. Its background part is

A=-(23/8+epsilon)N^(-3/2)
  -(2849/200-2epsilon)N^(1/2)
  +(9/8-epsilon)N^(5/2)+24+N^(-1/2)/200.

At epsilon=10^-6 its first five lapse derivatives at N=1 are

(801/100, 0, -749377/250000, 1077891/31250,
 -31494387/200000).

The linear and quadratic invariant polynomials, and the
two rational bounce mass functions, are stored explicitly
in center.py and the report. Each invariant coefficient
is compared with the original full Hamiltonian, not only
with its value at zero invariants. The generic S5 stationary
series identities also apply to this A+L+Q2 decomposition,
with A''=-2J_e; they retain the second-order lapse correction
in the fourth-degree reduced Hamiltonian. This remains
invariant degree, not a physical scattering amplitude.
