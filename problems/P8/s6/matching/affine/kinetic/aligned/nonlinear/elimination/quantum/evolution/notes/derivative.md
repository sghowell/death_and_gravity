# Uniform first time derivative of the finite vector integral

All quantities here are in tau=1 units until physical units are
restored. Retain exactly S6.50's fourth-order Gaussian preparation
and S6.53's matched prescription. Let I=[-1/2,1/2], m>=1000,
a in [1,Amax], Amax=25/16, omega^2=m^2+k^2/a^2, and
nu^2=m^2+k^2/Amax^2. Every estimate holds throughout I and for
every real k>=0. Endpoint derivatives mean continuous one-sided
limits. No state or frozen counterterm is changed.

## Actual differentiated reference bounds

The S6.50 reference is W=omega*S, S=1+P2/omega^2+P4/omega^4.
Its exact residual is rho_res=omega^-4 F. S6.50 already proves
S between 1/2 and 3/2, |W'/W|<4 and
|rho_res|<=C/omega^4, C=2,000,000.

The new calculation differentiates the actual rational P2,P4,
B2,B4,A,B using the complete fixed-comoving-momentum derivative,
including z'=-2Hz(1-z). Exact polynomial-box envelopes bound each
coefficient and its derivative on the whole I. Their denominators
are positive powers of 1+u^2. No time grid or truncated Hubble-jet
derivative is used.

Write F as a polynomial in these independent coefficients,
t=omega^-2 and S^-1, and write the S6.51 exact reference-tail
bracket the same way. For a polynomial sum c_alpha x^alpha,
its absolute-coefficient polynomial bounds its value. Differentiating
that positive polynomial against the stated absolute derivative
envelopes bounds the derivative by the chain and product rules.
The two polynomial recipes are independently replayed against the
frozen exact residual and subtraction identities.

Use |t'|<=4t and |(S^-1)'|<=4|S'|, with
|S'|<=|B2|t+|B4|t^2. For the physical readout coefficients,
|A|<=1, |B|<=3/2, |A'|,|B'|<2, |c1|<=2 and |c1'|<=5.
These last bounds follow directly from H=4u/(1+u^2), |H|<=2,
|H'|<=4, |z'|<=1, |alpha'|<=4/3, |beta'|<=28/27,
alpha<=4/9 and beta<=28/81. For example
|A_L'|<=4/3+4/9<2 and |B_T'|<=28/27+28/81<2.
The pressure weights have derivatives <=2/3.
The largest c1 derivative is at most 4+1=5.

The resulting residual derivative constants are below
86,409,927 and 91,121,821 for the two polarizations. Thus

    |rho_res'|<=C1/omega^4, C1=100,000,000.

The derivative of omega^-4 supplies 8|F| to this envelope.
For the reference-tail density R/(4a^3 omega^5), differentiation
gives [R'-(3H+5lambda)R]/(4a^3 omega^5), lambda=omega'/omega.
Here |3H+5lambda|<=16. The resulting bracket envelopes are below
20,376,082 and 19,921,102, so the common bound is

    |reference-tail density derivative|
        <=D1/(4a^3 omega^5), D1=100,000,000.

The report retains the rational envelopes, not only their ceilings.

## Remove the diagonal phase, then estimate mixing

Let f=(2W)^(-1/2) exp(-i theta), theta'=W, and write the exact
prepared solution v=Acal*f+Bcal*conj(f). The readout constraint
gives v'=Acal*f'+Bcal*conj(f'). The frozen CCR identity is
|Acal|^2-|Bcal|^2=1. Initially Acal=1, Bcal=0.

Set r=rho_res/(2W), eta'=r, eta(u0)=0, and
Acal=exp(-i eta)*a, Bcal=exp(i eta)*b. Exact variation of constants
then becomes

    a'=-i r exp(2i Phi) b,
    b'= i r exp(-2i Phi) a, Phi=theta+eta.

This is an exact change of variables, not deletion of the diagonal
term. Its harmless phase does not change either coefficient norm.
S6.50 gives |a|,|b|<=exp(J)<2, J=2C/nu^5. Also

    |r|<=C/omega^5, |r'|<=(C1+4C)/omega^5,
    Psi=Phi'=W+r>=omega/4, |Psi'|<=7omega.

The stated m floor verifies these inequalities explicitly. Therefore
g=r/(2Psi) has |g|<=2C/omega^6 and
|g'|<=(2C1+64C)/omega^6. Integrating the b equation once by parts,

    b(u)=-[g*a*exp(-2i Phi)]_u0^u
           +integral (g*a)' exp(-2i Phi) du.

Keep both endpoints. The length of I is one, and
|a'|<=|r| exp(J). Hence

    |b| <= exp(J)[2C1+68C+2C^2/nu^5]/nu^6
         < K/nu^6, K=4C1+138C=676,000,000.

The endpoint, bulk and feedback terms are all included. In particular,
2C^2/m^5<C and exp(J)<2. K/m^6<1.
The improvement is on actual mixing Bcal, not on the removable
diagonal phase. An undifferentiated bound on |Acal-1|+|Bcal|
would not justify the next step.

## Differentiated physical readouts

For either polarization let p=v'-d v. The exact canonical reduction
has U=d'+d^2, so v'=p+d v and p'=-d p-omega^2 v.
For a reference f the latter equation has the additional rho_res*f.
For any physical energy/pressure weights A,B define

    Q=[A|p|^2+B omega^2|v|^2]/(2a^3).

Its exact derivative has numerator

    [A'-(2d+3H)A]|p|^2
      +omega^2[B'+(2lambda+2d-3H)B]|v|^2
      +2omega^2(B-A) Re(p conj(v)).

The reference derivative has the additional numerator
2A*rho_res*Re(p_f conj(f)). Both identities are checked directly
using real and imaginary canonical coordinates, before integration.

The CCR and exact mixing expansion give differences of each squared
readout bounded by 6|Bcal| times the corresponding reference square;
the mixed real product difference is bounded by
6|Bcal| |p_f| |f|. Here |f|^2<=1/omega,
|p_f|^2<3omega, and |p_f f|<2. The diagonal phase cancels exactly.
Also |d|<=3, so the two differentiated diagonal weights are bounded
by 14 and 26, and |B-A|<=2. It follows that

    |Q'(v)-Q'(f)| <= [204|Bcal| omega+24|Bcal| omega^2
                        +2C/omega^4]/a^3
                    < 25K omega^2/nu^6+2C/omega^4.

The inequality 204/m+24<25 is checked. Summing three polarizations,
using omega<=Amax*nu and the exact integral
integral d^3k/(2pi)^3 nu^-4=Amax^3/(8pi*m), gives

    integral |exact/reference density derivative difference|
        < 40K/m.

The conservative integer 40 follows from pi>3 and
(75 Amax^2+6)Amax^3/24<40, with K>=C.
The reference-tail derivative integrates to at most D1/(72m^2),
as in S6.51. Thus for either physical energy or pressure,

    |d/du finite subtracted integral| <= 40K/m+D1/(72m^2).

These are uniform integrable derivative envelopes. Together with
the already established value integrals, they justify differentiation
under the integral and give a C1 finite state term on I. This does
not assert a second time derivative or all-order Hadamard regularity.
