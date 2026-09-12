# Full tensor reconstruction and the actual angular integrals

Choose phat=e3 only to perform the angular calculation. Use an orthonormal tracefree basis consisting of two tensor, two vector and one scalar component. Seven exact general-bilinear azimuthal identities, with independent coefficients for every detector and source basis tensor, reconstruct the full projector contractions. Off-diagonal spin components vanish under the azimuthal average, and each doublet has its correct common coefficient. This is checked before specializing the source and detector tensors.

The complete massive flat cubic and linear coefficients therefore combine into the invariant forms L, M and S in FORMULATION.md. The exact complete longitudinal moment is

integral_S2 |u|(n.D.n)(n.G.n) dOmega
=pi[2T+4V-W]/12.

Combining this with -3a P1_L/8 gives the actual curvature term

-p a(H_prime+2H^2)[2T+4V-W]/(1024pi^2).

It is determined by the original mode jets, not by choosing a finite curvature counterterm.

The complete centered conversion formula from S6.209 is

A1=(2pi)^-3 integral_S2 [
p^3 |u|(3-5u^2) f0/16 - p |u| f2c/2 ]dOmega.

Inserting the actual source-value and source-time-jet densities gives exactly all terms of FORMULATION.md. The original cubic coefficient is recovered at the same normalization and sign.

Independent finite-x angular checks use the actual four W8 frequencies and all ten physical readouts. For centered k=n+xP/2, define sigma=sqrt(1+pu x+p^2 x^2/4). The uncentered azimuthal geometry is evaluated at u_old=(u+px/2)/sigma and y_old=px/sigma. The 01 and10 factors acquire sigma^2, the11 factor sigma^4; the other four are scale-invariant. Eight literal physical azimuthal samples verify this exact transformation.

After coefficient extraction, the source-value second grade has angular degree at most six. Including the shell weight raises this at most to seven. Four Gauss nodes on[0,1] therefore integrate the even centered angular coefficients exactly; their moments through degree seven are separately checked in exact radicals. Full actual angular coefficients at twenty-four and thirty-two inverse-radius Cauchy resolutions agree.

The continuous extension at P0 follows because each direction form is uniformly bounded and every coefficient contains at least one factor p. This does not make the odd momentum magnitudes local or analytic at the origin.
