# Uniform squarefree energy derivative hierarchy

The normalized F=abc H3/W is homogeneous of degree2 in energies.
Set c=1 to construct it in QQ[a,b,r,v,u]. Ten dummy root exponents
are removed only after asserting they are zero and checking the exact
inverse representation map.

For each known denominator product D=product f_k^p_k, differentiation
uses the exact factored quotient rule. The common multiplier for a
variable x is L=product_(f_k,x nonzero) f_k and
 T=sum p_k*(partial_x f_k)*(L/f_k).
Thus the new numerator is L*(partial_x P)-T*P and the active powers
increase by1. Every L/f_k division has zero remainder.

After derivatives a,b, homogeneity gives degree0. Therefore at c=1,
 partial_c partial_a partial_b F
 =-a partial_a(partial_a partial_b F)
  -b partial_b(partial_a partial_b F).
The implemented Euler expression retains its exact common factors.
This proves the formula at all energies by homogeneity, not by sampling.

All denominator coefficients after division by
 S=1/(a+b)+1/(a+c)+1/(b+c) are positive. For every numerator monomial
the verifier either uses the same denominator monomial or locates two
supported exponent vectors symmetric about it. On positive variables,
 X^m <= (X^(m+d)+X^(m-d))/2.
It accumulates absolute rational coefficients and takes the maximum
coefficient ratio. There are no unsupported numerator monomials.
The sum of diagonal bounds plus twice off-diagonal bounds majorizes
the spatial Frobenius norm.

For lower derivatives the exact polynomial is homogenized to
QQ[a,b,c,r,v,u], with numerator degree denominator degree+2.
Setting c=1 is checked to recover the original polynomial exactly.
The same quotient rule bounds each first derivative divided by W
and each distinct second derivative divided by W*S. Homogeneous
midpoint shifts preserve energy degree.

All16 sector/basis cases have explicit executable caps:
 C1<=1000000, C2<=1000000, C3<721.
They are algebraic output inequalities, not chosen unknown physical
matching coefficients. The much looser C0=4800000000000000000000
comes from S318's complete n=3 weighted-current theorem.

An independent exact nilpotent ring in ea,eb,ec, modulo their squares,
retains the original vertices and original inverses. It supplies all48
constant/squarefree component coefficients for each reconstructed case.
Two additional configurations compare all96 such coefficients with
ordinary symbolic differentiation. Point checks calibrate the algebra;
positive-coefficient certificates prove the uniform inequalities.
