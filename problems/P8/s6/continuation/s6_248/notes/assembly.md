# The actual current operator: include the full heavy profile once

This is the existing CD-REG-AFFINE-ISO-QG2-H8A420 reference, on I=[-1/2,1/2], with a=(1+t²)², Proca mass1000, heavy mass squared n=10^200/512+2, and kappa=10^800. The full actual S240 SLE, S238/S239 parent and finite extensions, and S241 current classical coefficient functions are unchanged.

S241's coefficient sector already contains the entire heavy clock profile. S246's full clock response contains the pure heavy Gaussian response AND that profile. Consequently the correct scalar Euler assembly is

T_cl,QG2 + R_Proca,Gaussian + (R_S246,clock-P_H,clock)
 = T_cl,QG1 + R_Proca,Gaussian + R_S246,clock.

The equality uses T_cl,QG2-T_cl,QG1=P_H,clock at quadratic order. It is not a license to use the old QG1 coefficient comparison without checking the new coefficients. Adding the entire S246 response directly to the QG2 coefficient operator counts P_H,clock twice and has the nonzero defect P_H,clock.

To verify the equality, keep the full fixed functions rho_H(t), P_H(t) and set A_H=-P_H, B_H=-(rho_H+P_H)/2. For unnormalized heavy stress,

DeltaJ_H=(21delta²-3delta)A_H/2+(1-6delta)B_H,
T_H=(1+3delta)A_H-2B_H=rho_H-3delta P_H.

The complete mixed clock-profile Hessian for Q=2vI is

a³[2DeltaJ_H n_D n_G
 +3T_H(n_D v_G+n_G v_D)+9A_H v_D v_G].

It follows directly by taking the mixed derivative of

a³ N exp(3v) [1+2delta(N^-2-1)]^(-3/4)
 [A_H+B_H(N^-2-1)]

at N1,v0. This expression is used only for the exact degree-two jet: S240/S246 already verify that the ENTIRE fixed profile and ENTIRE parent R have precisely these value/first/second X jets. They are not replaced away from X1. Dividing the stress by kappa gives exactly the actual added S241 quadratic coefficient action. The code compares the actual full functions, not only formal names or numerical envelopes.

The full high-order parent source and every relevant clock/value jet through total degree two vanish at this reference. The S238 change and S239 finite extension have zero relevant clock jets. Hence the additional coherent heavy mean has no linear mixing with the four scalar/clock/M1 variables. Its independent massive KG block remains. This does not remove its nonzero Gaussian metric response or any interacting light/mixed loop.

Use adapted prepared coordinates, for q=P²/a²,

n=eta', zeta=w+H eta, v=zeta-delta eta',
b=c'+q eta, sigma=r+ell_matter eta.

Here b is the physical shift divergence, not the auxiliary chart variable of S221. The inverse is eta=I n, w=zeta-H eta, c=I(b-q eta), r=sigma-ell_matter eta, with the ORIGINAL lower zero germ.

The complete coefficient action is
-3v'²+(Jc+w_bg²/2-3Theta²)n²+6Theta nv'
+sigma'²/2+w_bg n sigma'-3ell_matter v' sigma
+2b(v'-Theta n)+ell_matter b sigma
+qv²+2E qnv-q sigma²/2+3Tcorr nv+9A v²/2,
with E=1-3delta and w_bg=-ell_matter E.
Jc,A,Tcorr are the WHOLE QG2 functions of S241. In particular DeltaJ_H is not added a second time to Jc.

Vary with the density-weighted adjoint W=D+3H after substituting all adapted variables. The full mixed row orders are(4,4,4,2), and their principal matrix is diag(-6delta²,0,0,-1). Every above-order coefficient vanishes. The full heavy profile has at most one derivative on either adapted leg, so it changes no displayed principal entry. Without the matter shift sigma=r+ell_matter eta, the matter row contains the nonzero coefficient ell_matter of eta''.

On I, delta>=32/125, so |(-6delta²)^-1|<=15625/6144. No division by H, Theta, E or external momentum is needed for this pivot. S241's explicitly rechecked current coefficient hypotheses supply the smooth finite-band classical comparison; neither the action's sign nor this pivot bound is a stability theorem.
