# Why the mixed functional has no extra reducible terms

Integrating the quadratic Dirac fields gives a bosonic action
S_B^(0)+h F, where F=-Tr log D_F and h counts loops. In the
bosonic loop expansion, the order-two term linear in this
order-one action is

    Gamma_2,F = (1/2) Tr[(S_B^(0)'')^-1 F''].

This is a statement about the effective action at fixed mean
background, not the connected vacuum functional at its
classical saddle. Those backgrounds differ. Leaving them
unconverted would retain spurious reducible terms.

The code independently computes a finite Euclidean Gaussian
example with S_0''=a>0, S_0'''=b, S_0''''=c, F'=f1 and F''=f2.
The exact second, fourth and sixth moments are 1/a, 3/a^2
and 15/a^3. The saddle calculation gives

    W2=-c/(8a^2)-f2/(2a)+5b^2/(24a^3)
        +b f1/(2a^2)+f1^2/(2a).

The classical-to-mean background shift is
d=(f1+b/(2a))/a. Its Legendre contribution yields

    Gamma2=-W2+a d^2/2
           =c/(8a^2)-b^2/(12a^3)+f2/(2a).

Both the two-fermion-loop dumbbell f1^2 and its mixed
fermion/bosonic-tadpole counterpart cancel. In arbitrary
finite indices the same operation cancels
(1/2)(F_i+t_i)D_ij(F_j+t_j), with
t_i=(1/2)C_ijk D_jk. Wick contractions and this bilinear
identity give the displayed trace with its factor one half.

Apply the finite-index derivation with a common regulator,
then interpret the indices as field labels and spacetime
arguments. The exact Gaussian H integral can be done first;
its determinant is Phi-independent and its full quartic
kernel remains. No fermionic quartic quantum vertex is
present in the original linear Dirac operator. At two
loops the remaining genuinely mixed topology is one
closed fermion loop with one contracted bosonic line,
dressed by background insertions.

This is a formal loop/diagram identity. The regulated
operator traces have not thereby been integrated, and
their errors or continuum counterterms are not bounded.
