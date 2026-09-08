# General-momentum selected-state mass bubble

Write p and r=K-p for the two comoving internal vectors and
c=p_hat dot r_hat. Directions are constant on the homogeneous
background, whereas their physical lengths and frequencies vary
with time. All v_T,v_L below are the unchanged exact S6.55 modes,
not their WKB references. Let p_L=v_L'-d_L v_L, with
d_L=H(1/2+q/omega^2). The physical orthonormal mode vectors are

    a^(3/2) U_T=(0,e_T v_T),
    a^(3/2) U_L=(sqrt(q) p_L/(m omega), i omega p_hat v_L/m).

The longitudinal temporal-spatial phase is essential. The exact
first-order equations v_L'=p_L+d_L v_L,
p_L'=-omega^2 v_L-d_L p_L imply the Proca divergence constraint.
At q=0, T and L obey the same equation and same selected preparation;
their sum tends to the isotropic three-spatial-component covariance.

Let W^+_{mu rho}(u,s,p)=sum_lambda U_lambda,mu(u,p)
conj(U_lambda,rho(s,p)). For V=m^2 diag(alpha,-beta,-beta,-beta),
J_density=a^3 V^{mu nu}W_mu W_nu/2. Wick's two connected contractions
give

    C^+(u,s,K)=a(u)^3 a(s)^3/2 integral d^3p/(2pi)^3
      V^{mu nu}(u)V^{rho sigma}(s)
      W^+_{mu rho}(u,s,p) W^+_{nu sigma}(u,s,K-p).

The two density volumes cancel the four a^(-3/2) mode factors.
Define pair amplitudes with these rescaled modes:

    A_TT=-m^2 beta (e_p dot e_r) v_Tp v_Tr,
    A_TL=-i m beta omega_r (e_p dot r_hat) v_Tp v_Lr,
    A_LT=-i m beta omega_p (p_hat dot e_r) v_Lp v_Tr,
    A_LL=alpha sqrt(q_p q_r) p_Lp p_Lr/(omega_p omega_r)
         +beta omega_p omega_r c v_Lp v_Lr.

Each connected integrand is sum A(u)conj(A(s))/2.
The transverse sums give angular weights 1+c^2 for TT and
1-c^2 for TL and LT; LL retains both terms and their interference.
The implementation checks all four against explicit polarization
vectors and against the four-index Wick contraction.

For K=0, c=-1, the mixed channels vanish and
A_LL=alpha z p_L^2-beta omega^2 v_L^2, as required by the
independent canonical insertion.

## Retarded sign and exact covariance comparison

Since H_int=-n J, the response is

    R_bubble(u,s,K)=i theta(u-s)[C^+(u,s,K)-C^-(u,s,K)]
                  =-2 theta(u-s) Im C^+(u,s,K),

where the last equality uses Hermiticity and this state's spatial
parity; C^- is the Fourier transform of the reversed current
ordering. Add the constraint contact
-delta(u-s) integral d^3p/(2pi)^3 alpha(u)^2 z(u,p)|p_L(u,p)|^2.
It is spatially local, so its momentum coefficient is independent
of external K. These formulas use coordinate-density outputs;
divide by a(u)^3 for the normalized physical output.

For x=(v,p), the symmetric covariance Sigma and symplectic matrix
J_2, an arbitrary real transfer F gives W=F(Sigma+i J_2/2).
With quadratic output A and source B, direct Wick contraction is
Tr(A W B W^T)/2. Its Kubo commutator equals
Tr(A F(delta M Sigma+Sigma delta M^T)F^T)/2,
where delta M=-J_2 B. The identity holds with arbitrary real
symmetric A,B,Sigma and arbitrary F; no stationarity is assumed.

## Flat frequency poles

With z=1-m^2/omega^2 and external Omega, the three-polarization
retarded flat kernel per momentum is

    [m^4 beta^2/omega+omega^3(beta+alpha z)^2/2]
      /[4omega^2-(Omega+i0)^2]-alpha^2 z omega/2.

The first numerator already includes the two transverse modes.
The common momentum measure is d^D p/(2pi)^D. For D=3-2epsilon,
the residues of the integrals of omega,1/omega,1/omega^3 are
respectively -m^4/32,-m^2/8,1/4 times 1/(pi^2 epsilon), obtained
independently from the radial Gamma integral. The coefficients of
Omega^0,Omega^2,Omega^4 in the complete response pole are

    3m^4(alpha^2+2alpha beta+5beta^2)/256,
    -m^2(alpha+beta)(5alpha+beta)/256,
    (alpha+beta)^2/512,

times 1/(pi^2 epsilon). All three equal the actual Lorentzian
timelike restriction F/(32pi^2 epsilon) of frozen S6.48, including
its quadratic-action to linear-response factor of two.
The omitted contact changes the residue by +3alpha^2 m^4/64.

These are pole checks only. Continuing the number of transverse
modes and the remaining D-dependent integrand affects the finite
part and must not be replaced by an early D=3 limit.
