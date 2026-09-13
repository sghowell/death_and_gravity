# Complete physical probes, remaining modes and the held-vector contact

Let every probe have support away from the common initial Cauchy
neighborhood. In unitary clock gauge choose

    N=1+epsilon eta(t),
    a_phys=a_ref exp(epsilon v_phys(t)),
    chi=chi_ref+epsilon psi(t),
    a_hat=R(u,N^-2)^1/4 a_phys.

The fixed Hclock(u), both fixed stress profiles, kappa, zeta and all masses
are not varied. Put r1=R_N|1, r2=R_NN|1. For the first and second
derivatives of log(a_hat) along this physical path,

    l1=v_phys+r1 eta/4,
    l2=(r2-r1^2)eta^2/4.

Therefore a_hat derivatives are a*l1 and a*(l1^2+l2), Hhat derivatives
are d_t l1 and d_t l2, and the fixed-comoving q derivatives are
-2q*l1 and q*(4l1^2-2l2). All time derivatives of eta, v_phys and psi
needed through this order are present. Treating q as fixed physical
momentum or treating Hclock as varied would change the result.

For each parent family f, the jet of its jth lapse derivative is

    (f_j, f_(j+1)*eta, f_(j+2)*eta^2).

The entries are derivatives, not Taylor coefficients. Jet multiplication
has second component x2*y0+2*x1*y1+x0*y2. Fractional powers use the full
second chain rule. Lapse derivatives0..4 suffice because J and the other
Gaussian coefficients contain at most two lapse derivatives already.
Coefficient time derivatives act on the fixed functions and on the
physical probe jets. The complete18-row map is stored for

    D,Z,J,K,Theta,w,c,Lnv,Vvv,Vvs,C,Y,Yv,r,rN,dH,q,a_hat.

For the entire canonical6x6 Hessian H(p), physical.vertices lists

    H_first=sum_p H_,p p1,
    H_second=sum_p H_,p p2+sum_(p,r) H_,pr p1 r1.

The symmetric double sum is implemented with all diagonal terms and
twice every distinct pair. EVERY matrix entry is retained. Substituting
each p,p1,p2 simultaneously from the explicit18-row map gives the
complete first and second directional physical vertices. Polarization
gives two distinct probe directions. The remaining variables have no
implicit missing variation. The exact reference specialization reproduces
the full S251 scalar matrix, full longitudinal Proca matrix and zero
spurious reference cross block. Nonzero physical first mixed vertices,
off-reference matter Ward derivatives and lapse second time derivatives
are explicit controls.

The other five modes have full physical canonical Hessians

    each TT: diag(2a_hat N C3 P^2, N/(a_hat^3 M)),
    each Proca transverse: diag(N(P^2/a_phys+a_phys/zeta),N/a_phys),
    H: diag(N(nH*a_phys^3+P^2*a_phys),N/a_phys^3).

Both copies of each polarization are included. All three derivative
orders are stored, not only the reference operators. The H source has
clock degree8 and hence no background/fluctuation derivative of total
degree at most4 here; its entire free determinant remains mandatory.

## Holding W0 fixed is different from aligning it

For a general held temporal background Wbar, define

    num0=Wbar-3r dH,
    num1=-A0-r(3vdot-b)-3rN dH n,
    num2=-b A_L/P-3r v b-rN n(3vdot-b)-3R_NN dH n^2/2.

The entire quadratic normal square is

    Z(num1^2+2num0 num2)/2+(Z_N n+3Zv)num0 num1
      +(Z_NN n^2/4+3Z_N n v/2+9Zv^2/4)num0^2.

The independent Fourier test expands N^x=b sin(Px)/P, v=v cos(Px),
Wi=A_L sin(Px) and W0=Wbar-A0 cos(Px), keeping shift transport and all
lapse/volume factors before averaging. The aligned restriction num0=0
recovers the full previous square. The1/P pairing is for P>0 with the
existing physical-shift infrared domain, not a literal homogeneous
zero-mode constraint formula.

At the reference, differentiate before eliminating the auxiliaries.
The envelope theorem gives the normal-vector Hamiltonian vertex

    h_W=rN*n*(3Theta*n-3c sigma/2-pv/2)
          +pv*A_L/(2P)+(Z_N*n+3v)*P*pA,
    n=[Theta(pv+3c sigma)-w ps-Lnv v]/(2J).

Its entire canonical Hessian and reference substitution are stored. In
the centered product reference the cross scalar/vector means vanish,
but the first scalar quadratic term need not have zero expectation.
Centering the free vector does not solve its interacting normal Gauss mean.

The physical aligned background has S_first=0 and

    S_second=6r1 eta*d_t(v_phys+r1 eta/4).

For any differentiable full finite Hamiltonian,
H_aligned_second=H_held_second+H_W*S_second. The complete held-W second
vertex is therefore the aligned second matrix MINUS the displayed whole
normal-vector canonical matrix times S_second, with its reference
coefficient substitutions. For an effective action the corresponding
identity is Gamma_aligned_second=Gamma_held_second+Gamma_W*S_second.
The finite Gaussian Hamiltonian contribution to Gamma_W is
-expectation(H_W); unconstructed nonlinear measure or other physical
counterfunctional contributions are not silently assumed zero.

These are background probes in a specified fixed affine-clock fluctuation
chart. A nonlinear change of the quantum fluctuation split also produces
the action-onepoint times the second field-map derivative. It requires its
own complete measure and contact bookkeeping. This checkpoint does not
identify the present Gaussian chart calculation with an arbitrary
covariant physical stress or a completed clock Ward current.
