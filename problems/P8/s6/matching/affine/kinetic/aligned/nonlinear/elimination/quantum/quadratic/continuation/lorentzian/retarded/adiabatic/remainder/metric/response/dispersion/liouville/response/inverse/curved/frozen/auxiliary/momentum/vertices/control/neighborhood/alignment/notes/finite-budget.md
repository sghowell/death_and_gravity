# Positive majorants and finite same-scale model-change budget

Use the unchanged seven free mode columns, original vector covariance,
two scalar charts, hard-transfer masks and fixed-total-momentum fibers
of S6.78. No nonlinear equality of canonical momenta is assumed. The
comparison uses the identity between the two specified canonical phase
charts, whose tangent free variables and quadratic Hamiltonian coincide.
It is a defined finite-order transition comparison in those variables,
not an invariant all-orders matching theorem.

The full S6.77 positive physical Hamiltonian majorant is monotone in
the common lapse-derivative bounds B0,...,B4. Let each B be multiplied
by R>=1. Write L_k=B_k*L and Q_k=B_k*Q, where L begins at physical
degree one and Q at degree two. The stationary majorant uses

n1=20*L1,
F2=B3*n1^2/2+L2*n1+Q1,
h=B0+L0+Q0+10*L1^2
  +B3*n1^3/6+L2*n1^2/2+Q1*n1
  +B4*n1^4/24+L3*n1^3/6+Q2*n1^2/2+10*F2^2.

The actual volume series multiplies h and the positive moving-boundary
majorant is added; neither depends on B. Using independent positive
coefficients for each physical degree in L,Q,volume and the moving
terms, exact expansion verifies positivity of every coefficient.
The degree-three coefficient has R degree at most four. The
degree-four coefficient has R degree at most six; the latter includes
the full F2^2 correction and is not just the direct quartic contact.
Thus h3_new<=R^4*h3_old and h4_new<=R^6*h4_old. The new twelve-row
inventory is bounded by the old fourteen-row majorant with the larger
B; zero rows do not invalidate a positive majorant.

Take the proved common R=10^9. All phase columns, spatial derivatives,
inverse transfers, fixed-center volume, Wick factors, Schur measures
and 84 connected exchange contractions are the previous ones. If C3
and C4 are the S6.78 scale-free norm numerators, then

C3_new <= R^4*C3,
C4_new <= R^8*C4.

The R^8 accommodates two cubic vertices and also dominates the R^6
quartic contact. By the triangle inequality, comparing both models at
the same scale L=M*tau gives model-change bounds

||Delta cubic|| <= (R^4+1)*C3/L,
||Delta connected quartic tree|| <= (R^8+1)*C4/L^2.

Name the new sufficient example **S6.80-P400**, L=10^400, with the
same fixed dimensionless vector mass 1000 and zeta=10^-6. Exact
rational checks give a cubic change below 10^-14 (approximately
3.426*10^-15) and a connected quartic change below 10^-150
(approximately 1.371*10^-153). Each new block is below 10^-3.
Decimals are descriptive only; the certificate uses rational bounds.
This enormous sufficient hierarchy is deliberately conservative and
not an optimum or a phenomenologically motivated scale proposal.

The observable is still only the finite-window cubic few-particle
blocks and hard-masked connected 2-to-2 tree block through H4 and
H3 squared. There is no claim about a full Fock operator norm,
vacuum/spectators, all interactions, forward transfer, loops, omitted
operators, thresholds, Wilsonian cutoff or UV matching. These are
not the former L=10^24 or L=10^353 numerical examples.
