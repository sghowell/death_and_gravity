# Complete finite-transfer canonical swap and outer reduction

Write J for the current positive Jc only in this note's algebra. The unchanged S220 phase Hamiltonian is
Hc=ps²/2-ell pv sigma/2+(q/2-3ell²/4)sigma²-qv²-9Av²/2+Lc²/(4J),
Lc=Theta pv-w ps+3Theta ell sigma-(2Eq+3Tcorr)v.
The weighted canonical equations are field'=H_momentum and momentum'+3H momentum=-H_field. The lapse is n=Lc/(2J). No Theta division occurs in this original system.

For the central chart set b=-pv/(2q). Integration of the full a³-weighted canonical term gives
pv vdot -> 2q v bdot+2Hqvb,
because q'=-2Hq. Dropping the second term changes the complete dynamics and its principal gradient.

Let x=(v,ps), y=(b,sigma), M=diag(2q,1), and Ctransport=diag(2Hq,0). Substitute pv=-2qb into the full Hc, and denote its x,x; x,y; y,y Hessians by XX,XY,YY. The canonical action is
x^T M ydot+x^T Ctransport y-Hc.
Set S=Ctransport-XY. On the certified high-q central domain XX is positive definite and
x=XX^-1(M ydot+S y),
Lcentral=(M ydot+S y)^T XX^-1(M ydot+S y)/2-y^T YY y/2.
Thus K=M^T XX^-1 M, B=M^T XX^-1 S, D=S^T XX^-1 S-YY. Full symbolic substitution and an independent stationary-action solve agree. No A,Tcorr or finite-q term is dropped.

Exactly,
detXX=[(2Eq+3Tcorr)²-(2q+9A)(2J+w²)]/(2J).
With z=1/q, Delta=(2E+3Tcorr z)²-(2z+9Az²)(2J+E²ell²), this is q²Delta/(2J). XX_ps,ps=1+w²/(2J)>1. The next note proves Delta>1/8 on the actual domain. The central chart uses no Theta denominator.

The actual high-q kinetic and gradient matrices are
K0=[[2J/E²+ell²,ell],[ell,1]],
G=[[2C/E²,ell],[ell,1]].
They reproduce twice the earlier core gamma module's half-Hessian convention. This old qualitative principal algebra is explicitly acknowledged as an input.

The outer chart substitutes n=(vdot+ell sigma/2)/Theta in the full S220 action. Its complete kinetic matrix is
K=[[ (2J+w²)/Theta²,w/Theta],[w/Theta,1]],
and its principal gradient is
G=[[2C/Theta²,w/Theta],[w/Theta,1]].
This chart is used only where |Theta|>=3/10. Both chart characteristic determinants give the same speeds.

For either exact L=ydot^T K ydot/2+ydot^T B y+y^T D y/2, retain the full time derivatives of all coefficients, including q'. Define Omega=B-B^T and R=B'+3HB-D-qG. The exact Euler operator is
K ydd+(K'+3HK+Omega)ydot+(qG+R)y.
Central K-K0 is O(1/q), while Omega and R are O(1). Outer K,Omega,R are independent of q. This is checked on complete rational expressions, not inferred from the principal pair.

The independent dynamical tests start from a separately entered smooth, nonzero-retuning Hamiltonian, differentiate the original phase coordinate map in time, and compare the complete transformed canonical generator and its derivative against both Euler operators. These algebra fixtures do not approximate the extremely small physical stress numerically. The physical coefficient bounds come from the exact fixed-profile proof.

The central chart really fails at some finite q: at the bare bounce detXX=0 at q=152/25. The original phase Hamiltonian remains finite there. Low momentum must therefore use that original regular system, not an unqualified extension of the b chart.
