# Leading integral, Laurent finite part and the exact correction

Put r=1/T, z=T/v and k=r z. The exact outer
Feynman parameter factor is

 J(k,epsilon)
  =[1-k^(1-epsilon)]/[(1-k)(1-epsilon)].

At mu=m, the regulated reference is

 F_D/(C/Q)=Gamma(epsilon) H(epsilon)
  integral_0^1 z^(2epsilon-1)
    (1-z)^(3/2-epsilon) A(k) J(k,epsilon) dz,

 H=e^(2gamma epsilon)4^(-epsilon)
       sqrt(pi)/(2 Gamma(3/2-epsilon)),
 A(k)=(1-k)^(-2).

Set r=0 only to define a separate reference
integral, not to discard the actual correction.
The beta integral then gives

 F_D,0/(C/Q)
 =e^(2gamma epsilon)4^(-epsilon) sqrt(pi)/2
  [(3/2-epsilon)/(1-epsilon)]
  Gamma(epsilon)Gamma(2epsilon)/
  Gamma(5/2+epsilon).

Multiplying by 2epsilon^2 gives an analytic
normalized function B with B(0)=1. Its
logarithmic derivatives are -7/3 and
5+pi^2/3. Hence

 F_D,0/(C/Q)
 =1/(2epsilon^2)-7/(6epsilon)
   +47/18+pi^2/12+O(epsilon).

The recurrence Gamma(x+1)=x Gamma(x),
differentiated twice logarithmically, gives
psi_1(x+1)=psi_1(x)-1/x^2.
Together with psi_1(1/2)=pi^2/2, two shifts
give psi_1(5/2)=pi^2/2-40/9. The code applies
these shifts explicitly, not a changed
special-function implementation. The half
value also follows from the absolutely
convergent trigamma series by splitting
the even/odd reciprocal squares.

For the exact difference define h=A-1 and
L(k)=-k log(k)/(1-k). Its compact integral
D(epsilon) has

 D0=integral w h(k) dz,
 D1=integral w{[2log z-log(1-z)+1]h(k)
                         -A(k)L(k)} dz,

where w=(1-z)^(3/2)/z. The constant term
of Gamma(epsilon)H(epsilon) after its pole
is 2-4log 2. Thus the exact finite difference is

 (F_MS-F_MS,0)/(C/Q)
 =integral_0^1 w{
   [3-4log 2+2log z-log(1-z)] h(rz)
   -A(rz)L(rz)} dz.

D0 is generally nonzero and gives a simple
pole D0/epsilon in the difference. MS removes
that pole; it is not treated as a vanishing
mass-ratio contribution. The finite formula
retains its product with the full prefactor.

The subtracted compact integral is
holomorphic for |epsilon|<1/8. At z=0,
A J-J(0,epsilon) is bounded by a constant
times z plus z^(1-Re epsilon), so the worst
combined power is integrable. At z=1 the
power is at least 11/8. Parameter derivatives
add only integrable logarithms. Normal
convergence therefore justifies extracting
D0,D1 before the endpoint integrations.
