# Complete retained canonical Hamiltonian and regular local branch

Use the source conventions in source.md. Normalize the metric trace
momentum by p=2 tr(pi_metric)/(3V), let G=partial_i pi_W^i/V, and set
a_trace=-M/3. Spatial norms in the following expression use gamma, with
momenta divided by V. The normal Hamiltonian per V has the full raw form

    Hraw=N[(pK-Ltrace)|Kstar-TG+Eother],
    Kstar=[p-B+U r(T-c)]/[2a_trace+U r^2],

where

    Eother=2|pi_TF/V|^2/M+(p_M1^2+p_H^2)/(2U)
        +Cchi(|grad M1|^2+|grad H|^2)/2
        +U(n h^2/2-jh)-C3 R3
        +|pi_W/V|^2/(2 zeta Cchi)+zeta M |Fij|^2/4
        +Cchi |Wi|^2/2.

No potential, curvature, gradient, electric, magnetic or heavy-source
term has been omitted. The full density adds Ni H_i, where on the primary
surface

    H_i=-2 gamma_ij D_k pi_metric^(jk)+pi_W^j Fij
        -Wi partial_j pi_W^j+p_M1 partial_i M1+p_H partial_i H.

Off that surface include the lapse and temporal scalar primary-momentum
transport terms. This is a classical spatial-gauge generator statement;
its quantum regulator and gauge-fixing determinant are not calculated.

Define Gamma=1-3r^2/(2R). The exact identity

    Gamma-1/4=3(R-1/2)(2-R)/(2R)

gives a strictly positive trace-velocity pivot on the original
1/2<R<6/5 domain. The five shear, three electric and two matter velocity
blocks are also invertible there. For the joint trace/temporal equations,

    Kjoint=(p-B-rG)/(2a_trace),
    Tstar=r Kjoint+c-G/U,
    Hred=N[(p-B-rG)^2/(4a_trace)-Fhat+G^2/(2U)-cG+Eother].

The full auxiliary Hessian D=d_(N,T)^2 Hraw obeys

    D_TT=-N U/Gamma,
    D_NT|Tstar=-D_TT partial_N Tstar,
    (D_NN-D_NT^2/D_TT)|Tstar=partial_N^2 Hred.

The derivative on the right is at fixed canonical momenta and all
retained fields and spatial jets. It acts on the entire source and every
coefficient. Do not differentiate after putting velocities on shell.

## Actual complete source datum, not a diagnostic root

In homogeneous normalization V=1 and pi_volume=3p, the complete reduced
Lagrangian is

    L=-3D Hhat^2+3B Hhat+Z(m1^2+mh^2)/2+Potential(N),
    Hcan=-(pi_volume-3B)^2/(12D)
         +(pi_M1^2+pi_H^2)/(2Z)-Potential(N).

Here D=M/N, Z=U/N; Potential includes Fhat, the full heavy potential and
source. With momenta substituted only after differentiating,

    Hcan_N=-L_N,
    Hcan_NN=-2J,
    J=Cnn+3Theta^2/D-Z_N^2(m1^2+mh^2)/(2Z),
    Cnn=L_NN/2, Theta=-Hhat D_N+B_N/2.

The report checks the literal full homogeneous action against S256.
At its actual initial datum u=0, a_hat=1, N=1+10^-6, Hhat=h=mh=0,
use the entire positive M1 constraint root. The full R_u(0,N)=0 and
I(0,1)=0 with I_N(0,N)=0 imply B=B_N=0 there. Hence G=0, p=B=0,
Tstar=0 and partial_N Tstar=0. The complete normal auxiliary equations
vanish, and the full lapse Schur pivot is -2J with the actual source-bound
enclosure J>151/100. Consequently det D=2N U J/Gamma>0 (V=1).
Restoring a general volume supplies the corresponding nonzero V factors.

The complete Hamiltonian is algebraic in N,T and smooth in the retained
canonical variables and their finite spatial jets. The ordinary implicit
function theorem at this datum therefore gives one local regular branch
for those parameters; it applies cellwise or to a fixed finite regulator.
For a compact finite set of cells take the intersection of the local
neighborhoods. No uniform continuum radius, arbitrary inhomogeneous
solution, nonlinear well-posedness result or global uniqueness is inferred.

