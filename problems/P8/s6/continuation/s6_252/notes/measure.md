# Complete Gaussian constraints and the configuration density

Start from the entire S220/S241 scalar coefficient action, with current Jc,A,Tcorr,ell,E,Theta and q. Let W=kappa*a^3 and use the physical momenta conjugate to(v,sigma), not their normalized p variables. Both auxiliary variables are retained: lapse n and nonzero-momentum shift divergence b.

The complete two-velocity Hessian is A_vel=W diag(-6,1). Its invertible Legendre transform gives auxiliary equations C_n=partial_n H and C_b=partial_b H, with

D_aux=partial_(n,b)^2 H=W diag(-2Jc,-2/3).

The second equation is C_b=(Pi_v-2Wb)/3; C_n includes Theta Pi_v-w Pi_sigma+W[3Theta ell sigma-(2Eq+3Tcorr)v]-2WJc n. All lower terms remain. The whole reduced Hamiltonian equals W times the S251 normalized Hamiltonian with its momenta divided by W.

On the extended canonical phase space, use constraints(p_n,p_b,C_n,C_b). Their full bracket is M=[[0,-D_aux],[D_aux^T,K]], where K_12=-W(2Eq+3Tcorr)/3 is generally NONZERO. Direct block multiplication gives inverse[[D_aux^-T K D_aux^-1,D_aux^-T],[-D_aux^-1,0]] and det M=(det D_aux)^2. Thus the positive constraint density is |det D_aux|=4W^2 Jc/3. It exactly cancels the auxiliary delta-function integration Jacobian. The remaining(v,sigma,Pi_v,Pi_sigma) bracket is canonical, because the lower-right block of M^-1 is zero and physical functions commute with the primary auxiliary momenta.

This is a finite-regulator phase-space identity for the already spatial-gauge-fixed QUADRATIC reference. It is not a proof of a full nonlinear covariant/BRST measure or a zero contribution from every parent gauge/connection determinant. The original coefficient action need not be stationary by itself. No nonlinear saddle or physical quantum constraint is substituted for this specified Gaussian system.

## Independent configuration integration

The auxiliary Hessian of the complete Lagrangian has determinant -4W^2 Theta^2. On an outer chart with Theta nonzero, exact configuration elimination gives det K_red=2W^2 Jc/Theta^2. The block determinant identity, including all velocity/auxiliary mixing, is

det K_red=det A_vel * det D_aux / det L_aux.

Consequently the full unreduced Gaussian configuration density, before n,b integration, is W^2 sqrt(8Jc), up to a common coefficient-independent normalization and consistently continued Fresnel phase. Multiplying only by sqrt(|det A_vel|) would miss a Jc-dependent factor. After the auxiliary Gaussian integral the density is W sqrt(2Jc)/|Theta|, exactly the reduced canonical momentum integral. The determinant cancellation must precede any continuum regularization; separately regularized functional determinants need not obey an unqualified product rule.

For canonical fields sqrt(kappa)(v,sigma,n,b), the four-coordinate Jacobian changes the squared unreduced density to8a^12 Jc. For the two reduced canonical fields the outer squared density becomes2a^6 Jc/Theta^2. These are both the same physical normalization as S251.

At Theta0 this outer configuration integration is not used. The original phase-space auxiliary D_aux is unchanged there and Jc>1/100. The certified central canonical chart has det(W K_central)=8W^2 Jc/Delta with Delta bounded away from zero in its stated q>=4096 range. Low momentum remains in the original regular phase space. No literal P0 shift constraint is assigned by its Fourier extension.

## Charts and boundaries

For a time-dependent canonical map y=Tz, the symmetric phase action -z^T Omega z'/2-z^T H z/2 gives H_new=T^-T H T^-1-Omega T'T^-1. Its ordinary p dq convention additionally has boundary(p_old q_old-p_new q_new)/2. The central swap is a Fourier-type canonical map, not a position-only relabeling; the symmetric-boundary shear has its full quadratic generating function. Transport the SAME initial Gaussian density and boundary phases on both CTP branches. A unit determinant alone does not authorize an instantaneous vacuum reset.

The finite-measure framework is consistent with [Bellorin and Droguett, section5.1](https://arxiv.org/pdf/1912.06749). Their Horava theory is not our action; every displayed present block is derived independently. [Langlois and Noui, sectionVI](https://arxiv.org/pdf/1512.06820) explain why a regular unitary chart must be checked. No statement here extends that chart through a vanishing clock gradient.
