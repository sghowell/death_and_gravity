# Ordered density Ward identities and complete chart correction

Use coordinate measure and define contravariant weight-one E by J(D)=integral E:hD, hD=Dg(D). The existing ADM current convention gives E00=A=-a^3 rho/2, Eij=B deltaij=-a P deltaij/2 and E0i=0.

The density Lie derivative is Lxi E=xi.partial E-(partial xi)E-E(partial xi)^T+(div xi)E. Its components are

    (Lxi E)00=eta A'-A eta'+A div chi,
    (Lxi E)0i=-B partial_i eta-A chi_i',
    (Lxi E)ij=eta B' deltaij-B symgrad(chi)ij
                         +B(eta'+div chi)deltaij.

The actual current is not zero. Conservation gives rho'+3H(rho+P)=0.

Let q(D,G)=D^2g(D,G) include every S214 spatial, lapse and shift chart term. When the covariant distributional variation is defined, its ADM response is R(D,G)=integral hD:deltaE[hG]+E:q(D,G). This chain rule does not make a retarded kernel a symmetric single-branch effective-action Hessian.

For an initial-identity gauge source, covariance of the Proca equation, Hadamard prescription and initial-state pullback gives deltaE[Lxi_G g]=Lxi_G E. Hence

    R(D,Gxi)=integral hD:Lxi_G E+E:q(D,Gxi).

Independently differentiate conservation with a fixed detector vector xi_D. The differentiated time flux is
2 xi_D^alpha[deltaE^{0nu}g_alpha nu+E^{0nu}hG_alpha nu].
It vanishes at the right because xi_D=0 and at the left because hG and the prepared current tangent vanish. Therefore

    R(Dxi,G)=integral -E:Lxi_D hG+E:q(Dxi,G).

Bilinearity yields the ordered reconstruction

    R(D,G)=R(Dsyn,Gsyn)+R(D,Gxi_G)+R(Dxi_D,Gsyn).

Do not interchange the last two terms or erase q. These identities concern the original covariant renormalized current, not an arbitrary-diffeomorphism Ward identity for a sharp computational momentum band. They do not prove the differentiability, finite matching or norms of the missing scalar kernels; reconstruction remains conditional on those inputs.

The full parent source is not deleted. S176's S=DS=0 identity and zero Gaussian mean of the linear second-source insertion supply only the conditional first metric-response bridge to ordinary connected Proca at this clock. Clock/scalar dynamics and any fixed scalar retuning must be varied separately.
