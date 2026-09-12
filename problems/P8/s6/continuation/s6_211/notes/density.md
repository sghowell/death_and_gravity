# Complete spatial difference and angular reconstruction

The seven exact azimuthal contractions are the original00,01,10,11,TL,LT,LL projector forms. No common inverse phase is substituted for the four sector phases. The existing general tensor reconstruction is replayed as a polynomial identity with independent five-component detector and source tensors. It proves diagonal tensor/tensor, vector/vector and scalar structure for every geometric contraction, not merely for one chosen polarization.

The unshifted azimuthal geometry uses y=px and u=n.e3. Multiply its exact degree4 expansion by each actual source-derivative amplitude row before integrating. The azimuth is already averaged; the remaining physical Fourier measure is

    integral_-1^1 du/(4pi^2).

Every coefficient is polynomial in u, so its integral is the exact sum of its even monomial moments2/(degree+1). Each inverse-radius degree d has angular degree at most d+4; five Gauss-Legendre nodes integrate the maximum degree8 exactly in the independent angular check.

Subtract only P0 from the fully summed angular coefficient. Among the105 channel/source/degree entries, tensor has4 nonzero differences, vector3, scalar4: eleven total. The vector source-value log difference vanishes; this is why counting four per channel would be wrong.

The only nonzero slots are j0/r0/d2 in all channels; j0/r0/d4 in tensor/scalar; and j2/r1/d2,j2/r2/d2 in all channels. Their invariant formulas appear in FORMULATION.md. All other differences vanish by the exact coefficient computation, including the j1/d3 difference. That endpoint itself has a nonzero homogeneous coefficient and remains in the full current.

At a fixed P the one-ball radial measure makes total grade q contribute r^(3-q)dr. The q0 difference is zero and all odd grades integrate to zero. Direct integration therefore gives

    DeltaU_ball=DeltaF2*(K^2-m^2)/2+DeltaF4*log(K/m).

The derivative with respect to K is exactly K DeltaF2+DeltaF4/K and the result is zero at K=m. Thus the finite lower term-m^2 DeltaF2/2 is retained, not hidden by asymptotic notation. The original two-leg sharp band is restored only through the unchanged complete S209/S210 conversion; this identity does not replace that regulator.
