# Complete canonical state and immutable reference prescription

This checkpoint concerns the local quadratic coefficient sector of the SAME QG2-H8A420 action. The full current S241 functions Jc,A,Tcorr, the original M1 charge and the actual scale factor remain. The classical first variation alone need not vanish: the retained Proca/H mean terms complete its existing stationarity. The new free state is a reference for additional quantum calculations, not a claim that those additional mean terms have vanished or that the already computed nonlocal response has been quantized. Formal loop orders must retain the known counterterm grades.

## Full physical momenta and both charts

Let Z=(v,sigma,pv,ps) be the original normalized variables, Omega=[[0,I],[-I,0]], and Dp=diag(0,0,1,1). The whole normalized Hamiltonian h is the S241 expression, with EVERY current coefficient substituted by its complete fixed function. Its equation is Z'=[Omega Hess h-3H Dp]Z. The physical momenta are kappa a^3(pv,ps), so the normalized equal-time bracket is Omega/(kappa a^3), not Omega.

Set z=C Z, C=sqrt(kappa)diag(I,a^3 I). The physical canonical pair is (sqrt(kappa)(v,sigma),sqrt(kappa)a^3(pv,ps)). The generator A_c=C A_Z C^-1+C'C^-1 equals Omega times the Hessian of kappa a^3 h(C^-1 z). The volume derivative cancels exactly the weighted momentum damping. Hence its exact full fundamental matrix S(t,s;P) is real symplectic for every finite momentum, including the continuous extension at zero. The full kappa normalization cancels from A_c, but is restored at BOTH field endpoints.

In the outer chart Q=sqrt(kappa)(v,sigma) and Pi=a^3(K Q'+B Q). Its exact canonical momenta coincide with those above after the full lapse/shift elimination. In the central chart Qb=sqrt(kappa)b=-Pi_v/(2a^3 q), Pib=2a^3 q Qv. The other pair is unchanged. This is an exact symplectic transformation with tau=2a^3q and tau'=H tau. Its time connection is indispensable. The central chart is used ONLY in its certified q>=4096 range; the original Hamiltonian supplies every low momentum.

For either complete reduced L=a^3[Q'^T K Q'/2+Q'^T B Q+Q^T D Q/2], write B=S_B+A_B into symmetric and antisymmetric parts. Subtracting the exact time boundary (a^3 Q^T S_B Q/2)' yields canonical momentum p=Pi-a^3 S_B Q and retains the gyro A_B. The complete symmetric potential is qG+V, V=S_B'+3H S_B-D-qG. The resulting Hamiltonian is

H_clean=(p-a^3 A_B Q)^T K^-1(p-a^3 A_B Q)/(2a^3)+a^3 Q^T(qG+V)Q/2.

Every q'=-2Hq and coefficient derivative remains. Exact four-by-four algebra checks the whole canonical generator against this expression in both charts.

After BOTH symmetric boundaries are retained, the central-from-outer transformation is T_clean=S_c T_swap S_o^-1. In the high-momentum scaled phase diag(sqrt(k)I,k^-1/2 I), it and its inverse are classical order-zero matrix symbols with leading map diag(-E/Theta,1,-Theta/E,1). The overlaps have nonzero E and Theta. Thus there is no singular state change at either switch +/-3/16. Uncleaned momenta can contain order-q boundary shears; counting their powers as physical ultraviolet excitation would be wrong.

## Specified coupled Gaussian selection

Fix t*=-1/2 and f(t)=exp[-1/(1-1024(t+7/16)^2)] for |t+7/16|<1/32, zero otherwise. This smooth nonzero bump lies strictly in the certified outer region. No arbitrary unspecified time or per-mode oscillator basis remains.

For a real classical canonical trajectory of the full physical scalar system define

e_t(z)=a^3[Q'^T K Q'+Q^T(qG+K)Q]/2.

The extra positive K term is a CHOSEN PREPARATION NORM at mass scale1. It is not the stress tensor, a mass correction to the physical action, or a retuning of the candidate. Write e_t=z^T E_t z/2, where Q' is obtained from the exact canonical generator, and set

M(P)=integral f(t)^2 S(t,t*;P)^T E_t(P) S(t,t*;P) dt.

For every finite P, E_t is strictly positive on the support: K,G>0, q>=0, and the exact Legendre map is invertible there. The continuous invertible S and nonzero bump give M>0. Smooth dependence of the exact ODE on P^2 gives a smooth positive M at P0 as well. This is not a positive-Hamiltonian or no-growth assertion on the whole slab.

Define p=sqrt(det M)>0, s=sqrt(-Tr[(Omega M)^2]/2+2p)>0 and

V(P)=[p M^-1-Omega M Omega]/(2s).

This formula preserves ALL off-diagonal mixing. In two modes it equals one-half M^-1/2 sqrt[-(M^1/2 Omega M^1/2)^2] M^-1/2. To prove it, Williamson-conjugate M to diag(nu1,nu2,nu1,nu2), nu_i>0. Both expressions become I/2 and transform covariantly back. Equivalently, the Hamiltonian matrix satisfies A^4+(nu1^2+nu2^2)A^2+nu1^2 nu2^2 I=0. The formula has no division by nu1-nu2 and remains regular at preparation-frequency degeneracies.

It follows that V+iOmega/2>=0 on ALL complex vectors and V Omega V=Omega/4. In a Williamson basis the energy is a sum nu_i(a_i^dagger a_i+1/2). Its unique zero-mean minimum is the pure joint vacuum, of energy(nu1+nu2)/2. This proves positivity, purity and minimality, not merely a real quadratic inequality. A real positive covariance smaller than I/2 fails the complex uncertainty condition and is an explicit rejected example.

The whole unequal-time canonical two-point matrix is

W_z(t,s;P)=S(t,t*;P)[V(P)+iOmega/2]S(s,t*;P)^T.

It is a bisolution, has the exact commutator, and is positive on arbitrary complex compact test vectors by integrating their pullbacks against V+iOmega/2. Recover the original fields and momenta with C^-1 at their OWN endpoint. The determinant/profile functions are not reselected on varied histories.

## Two tensor modes and direct product

The existing full current classical TT action has gamma=2h/sqrt(kappa) and each Frobenius-normalized physical polarization has action a^3(h'^2-qh^2)/2. Its canonical pair is(h,a^3h'), generator[[0,a^-3],[-a^3q,0]]. Apply the SAME f to positive selection form a^3[h'^2+(q+1)h^2]/2. Its complete two-by-two Gramian M_T is positive even at zero momentum, and the unique pure covariance is sqrt(det M_T) M_T^-1/2.

Both tensor polarizations use that covariance. The spatial covariance is obtained with the source-pinned coordinate-free TT orthogonal projector, whose norm is1. It is independent of the polarization frame chosen locally; no nonexistent globally smooth frame is asserted. The exact scalar/tensor classical sectors decouple, so their product gives four physical free modes. A coincident tensor/scalar speed1 creates no mixing between these exact rotation sectors. No vector gauge degree of freedom is added as a physical oscillator.

The original H and Proca Gaussian states and classical M1 mean remain. Tensor finite local/nonlocal quantum response is NOT deleted: it belongs to the omitted contributions beyond this chosen free-reference construction. A literal homogeneous lapse/shift constraint sector is not specified by a measure-zero Fourier extension.

## Holding the state fixed

The bump is a reference-state selection prescription, not a claimed operational measurement protocol. It defines Cauchy data at t* using the fixed reference coefficients. Subsequent perturbations keep these data fixed, with no new minimization, future boundary condition or state reset at a chart transition. This distinction also applies although the selection interval lies after t* on the fixed reference. Any interacting-state, physical loop, or live off-reference generalization must prove its own hypotheses.
