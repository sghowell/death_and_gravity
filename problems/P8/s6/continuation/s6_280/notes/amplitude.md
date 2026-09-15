# Independent literal Einstein and matter derivation

Use g=eta+2h/sqrt(kappa) and the original -kappa sqrt(-g)R/2.
After its exact divergence the density is
-kappa sqrt(-g)g^(mn)(Gamma^r_ml Gamma^l_nr-Gamma^r_mn Gamma^l_rl)/2.
For the three-plane-wave coefficient the total momentum is zero, so
the divergence vanishes exactly; no cosmological canonical boundary
is being dropped here.

At kappa1, g=eta+2 sum_i z_i H_i. Write Gamma=Gamma1+Gamma2 and
sqrt(-g)g^-1=eta+sum_i z_i[tr(eta H_i)eta-2eta H_i eta]+...
The code includes every Gamma1/Gamma2 and density/Gamma1/Gamma1
permutation containing each of the three distinct legs once.
An independent transverse-traceless quadratic wave checks the canonical
kinetic sign and normalization before this cubic vertex is used.

Put the outgoing physical gravitons along opposite z directions, with
incoming negatives k1=(-e,0,0,-e),k2=(-e,0,0,e). The incoming scalar
momenta are a=(e,p sin(theta),0,p cos(theta)),b=(e,-p sin(theta),0,-p cos(theta)),
m^2=e^2-p^2. Contract the cubic with the entire conserved scalar stress
through its de Donder projector. Add both scalar exchange diagrams and
the Phi^2 h^2 contact obtained from the complete quadratic expansion
of sqrt(-g)[g^(mn)partial_m Phi partial_n Phi-m^2 Phi^2]/2.

With D=e^2-p^2x^2, all four physical helicities give

    A_same = m^4/(kappa D),
    A_opposite = p^4(1-x^2)^2/(kappa D).

Each helicity tensor is epsilon epsilon with null transverse epsilon
and epsilon.epsilon*=-1. The full16-by16 polarization sum is checked
against the rank2 TT projector. No tensor product containing extra
scalar or antisymmetric polarization states is sewn.

A separate covariant factorization check uses t1=2a.k1,t2=2a.k2,s=2k1.k2
and numerator
N=2(epsilon1.b)(epsilon2.a)t1+2(epsilon1.a)(epsilon2.b)t2
  +(epsilon1.epsilon2)t1t2.
The literal result equals N^2/(kappa t1 s t2). Both Ward replacements
N(epsilon1=k1)=0 and N(epsilon2=k2)=0 are checked separately.

This is consistent with the covariant KLT organization in
[Bjerrum-Bohr et al.,section4](https://arxiv.org/html/1908.09755)
and the Compton factorization in
[Bjerrum-Bohr et al.,appendixA.2](https://arxiv.org/html/1609.07477).
Their printed spinor prefactors are not adopted as a normalization
axiom: here the full action, propagator and vertices fix it independently.

At threshold the two same-helicity amplitudes are m^2/kappa and the
opposite pair vanish. The m=0 formula is only a diagnostic limit; the
original field retains mass1.
