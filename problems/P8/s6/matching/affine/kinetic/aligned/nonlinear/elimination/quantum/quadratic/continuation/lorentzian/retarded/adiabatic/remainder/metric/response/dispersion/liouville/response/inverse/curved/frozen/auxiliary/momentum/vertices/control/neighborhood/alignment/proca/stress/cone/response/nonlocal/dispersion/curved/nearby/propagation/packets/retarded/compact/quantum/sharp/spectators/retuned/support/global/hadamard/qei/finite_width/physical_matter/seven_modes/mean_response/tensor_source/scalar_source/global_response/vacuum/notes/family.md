# Actual reduced classical family and the smooth branch

Use the repository +--- convention X=(d phi)^2, u=phi/tau,
h=(1+u^2)^3 and dimensionless action normalization
kappa=M^2 tau^2. The retuned tree scalar function is the full
p8_affine.dictionary scalar_F with x=-X plus
(X-1)^2/(200h^2). It is quadratic in X with denominator
200(1+u^2)^12; its actual constant at (u,X)=(0,0) is -28/25.
No quantum tadpole/profile function is included or reset.

Let r(z)=0 for z<=0 and exp(-1/z) for z>0, and

    B(X)=r(X-1/4)/(r(X-1/4)+r(3/4-X)),
    R=1+B(X)(X-1)/h.

Define F2=-R/2, A1=A2=K=0 and recompute the Ia completion:

    A3=R_X/X,
    A4=-R_X/X - 7 R_X^2/(4R),
    A5=R_X^2/(RX).

The last denominator means R times X. In the vacuum region
X<=1/4 define A3=A4=A5=0 directly and F2=-1/2. The C-infinity
flatness of the bump joins these formulas. This is not the
undefined product of zero with the old singular coefficient.

For X>=3/4 all full classical coefficient functions equal
the old retuned functions, for every u. Equality includes
the lower-order scalar function, not only the principal symbol.
The physical metric, free chi and ordinary Proca couplings
are unchanged. At kappa=1 the vacuum scalar kernel is

    Fv=Z X/2 - Z m^2(u-u0)^2/2 + lambda Z^2 X^2.

For general kappa multiply the last coefficient by kappa.
Then kappa Fv in Phi_bar=sqrt(kappa Z)(u-u0),
Y_bar=kappa Z X is exactly
Y_bar/2-m^2 Phi_bar^2/2+lambda Y_bar^2. Thus a fixed-mass,
fixed-interaction decoupling family is explicitly specified,
rather than silently holding the wrong bare coefficient fixed.
The switched lower kernel is B Ftree+(1-B)Fv.

For all real u,X, R>=83/164. If 1/4<=X<=1/2, use B<=1/2
and 1-X<=3/4 to get R>=5/8. On 1/2<=X<=3/5, the switch odds
at 3/5 are exp(80/21)<exp(4)<81, since e<3, so B<81/82
and R>83/164. For X>=3/5 use R>=min(1,X)>=3/5.
The vacuum plateau has R=1. Thus both tensor coefficients
equal positive R and |F2|>=83/328. This does not prove the
scalar kinetic determinant or a propagating heavy gap.

The necessary exceptional matter relation is
X A3+2F2_X=0. It holds identically in the completed family.
The tempting separate taper A3=B/(hX) misses
(X-1)B'/(hX). Substitution in the actual regular free-M1
source matrix gives

    f-g=-2X(X-1)B'/(h Theta),
    det(K-G)=-Y X^2(X-1)^2 B'^2/(h^2 Theta^2)<0

on a regular rolling transition configuration with Y>0 and
Theta nonzero. At u=0,X=1/2, B=1/2,B'=8,R=3/4,R_X=-7/2:
the correct A3 is -7, the naive A3 is 1, and the naive
difference determinant is -4Y/Theta^2. These are conditional
matrix facts, not an existence theorem for such a background,
a ghost proof or a universal UV no-go.
