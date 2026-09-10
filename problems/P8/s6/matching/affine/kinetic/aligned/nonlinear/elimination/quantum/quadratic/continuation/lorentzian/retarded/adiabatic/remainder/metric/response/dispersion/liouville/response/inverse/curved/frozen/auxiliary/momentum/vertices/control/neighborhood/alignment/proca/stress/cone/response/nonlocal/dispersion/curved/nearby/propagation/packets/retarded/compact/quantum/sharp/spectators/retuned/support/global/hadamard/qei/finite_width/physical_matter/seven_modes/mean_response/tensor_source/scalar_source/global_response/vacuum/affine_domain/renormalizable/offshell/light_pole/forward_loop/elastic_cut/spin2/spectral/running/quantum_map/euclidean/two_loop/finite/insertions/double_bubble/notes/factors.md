# Finite products with complete heavy triangles

At one external channel define C=-L+g/(M-z), where
L=lambda4 and z is s, t or u. Let T_A and T_B be the
two integrated internal-heavy choices at the respective
outer vertices, including g and the loop measure.

The 32 factorizing raw refinements are

    C (C I+T_A)(C I+T_B)/4.

Remove C T_A T_B/4, the eight finite refinements
already counted in S6.121. The selected raw group is

    [C^3 I^2 + C^2 I(T_A+T_B)]/4.

Its actual forest operations replace each remaining
divergent bubble factor by I_R=I-I0bar. Thus the
renormalized selected group is exactly

    [C^3 I_R^2 + C^2 I_R(T_A+T_B)]/4.

This identity holds before removal of the common
regulator. The finite limits follow from the
majorants below. The two independent loop integrations
are not confused with a single integral of a squared
vertex. Equality between T_A and T_B is not needed.

On the common first sheet,

    I_R(z)=-(16pi^2)^-1 integral_0^1
                  ln[1-x(1-x)z] dx.

For |z|<=3, the absolute logarithmic series is bounded
by ln(4)<2. The latter bound follows from the positive
partial sum 1+2+2=5 of exp(2). Therefore
|I_R|<2/(16pi^2). At the t-channel z=0 this factor
vanishes exactly, although the bound conservatively
retains all three channels.

For each finite triangle use the S6.122 explicit
Feynman shift with real radial y=q^2 and light
denominator gap delta=1/4. Every shifted heavy vector
has Hermitian norm below three. The square identity
gives

    Re[(q+w)^2+M] >= y/2+M-9 > (y+M)/2

at M>32. Hence each sum of two full heavy triangles obeys

    |T_A|, |T_B| <= 4g/(16pi^2)
       integral_0^infinity y/[(y+delta)^2(y+M)] dy.

The factor four includes two exchanges and the factor
two in the heavy-denominator bound. The anchored primitive

    M[ln(y+delta)-ln(y+M)]/(M-delta)^2
        +delta/[(M-delta)(y+delta)]

has the checked derivative and vanishes at infinity.
The exact integral is

    J=M ln(M/delta)/(M-delta)^2 - 1/(M-delta).

For delta=1/4 and M>=32, M^2/(M-delta)^2<2, giving
J<2ln(4M)/M. The discarded second term is negative.
These are all-radius absolutely convergent bounds,
not an expansion of the heavy propagator.

The same strict parameter and routing margins as the
parents justify the first-sheet continuation and a
neighborhood of the closed forward disc. Products
of these integrated holomorphic factors therefore
satisfy the ordinary Cauchy coefficient estimate.
