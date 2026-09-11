# Scalar source from both the covariance and fixed vacuum derivative

Introduce a light squared mass a as a background marker,
with a=1+GH near H=0. It is not a new moving OS condition.
Write T=T_D(1), I2=integral D1^2=-partial_a T_D(a)|1,
B0=B_HL,D(-1), r=-g partial_(pE^2)B_HL,D(-1), and

    S_a=partial_a S_D(a,a,M)|_(a=1).

Both light lines carry a, so
integral D1^2 B_HL,D(pE^2)=-S_a/2.
The first OS scalar covariance insertion integrates to

    Tinsert_scalar,D=-g S_a/2-g B0 I2+r T.

All three terms retain the common regulator.
The last is a correlated physical mass/residue reference,
not a discarded wavefunction insertion.

An independent derivation starts with the background
vacuum graphs L T_D(a)^2/8-g S_D(a,a,M)/4.
The reference counterterms stay fixed at a=1:

    deltaZ=-r,
    deltaMass=-L T/2+g B0-r.

Their determinant insertion at variable a is
(deltaMass-a deltaZ)T_D(a)/2 after its scaleless
constant trace is removed. Differentiating the sum
with respect to H, equivalently G partial_a at a=1,
gives G Tinsert_scalar,D/2. The raw local L derivative
cancels its physical mass-reference derivative exactly.

It would be incorrect to differentiate deltaMass or
deltaZ as if they were retuned at every H background.
That changes the specified theory and its one-point
function. A dedicated negative control detects this.

Adding the full fermion covariance insertion from
S6.148 and the proper cubic coefficient gives

    J2,H=Fin[-G(Tinsert_scalar,D+Tinsert_fermion,D)/2
             -deltaG1,D T_D(1)/2].

Here H denotes the fixed physical-Phi hybrid reference
when used as a scheme subscript, not a new field.
The overall source pole is removed after all proper
references are retained. Its finite value is fixed
by the original one-point-zero condition.
