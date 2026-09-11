# Full fixed profile and the actual conditional family limit

The separately named new action fixes
Deltaf(Phi,Y)=kappa0 DeltaF(Phi/sqrt(kappa0),Y/kappa0)
once, with the full coefficient

    DeltaF=-p_v+T(X)V(u,X),
    V=-p(u)+p_v-[r(u)+p(u)](X-1)/2,
    T=X^N/[X^N+(1-X)^N], N=1024.

Here r,p are the already source-defined reference functions,
not adjustable functions during variation. On the present
real canonical unit-jet class, |u|<=1/sqrt(kappa0) and
|X|<=4/kappa0. The S6.176 bounds imply
|V|,|V_u|<4epsilon and |V_X|<epsilon, epsilon=10^-770.
They apply because this u interval is inside its certified
compact reference slab. No complex analyticity of r or p
is assumed or needed.

The scalar constant -kappa0 p_v cancels the same-prescription
positive-frequency Minkowski vector vacuum constant. It
must be retained before discussing the gravitational
family. After that explicit cancellation, the remaining
coefficient is kappa0 T V, not zero. Its first possible
nonconstant vacuum field degree is2048, as established by
the actual full coefficient in S6.182.

For real |X|<=4/kappa0, even N gives a positive denominator.
Bernoulli's inequality yields
(1-4/kappa0)^(N+1)>1/2. Hence

    |T|<=2|X|^N, |T_X|<=2N|X|^(N-1).

The exact derivative is
N X^(N-1)(1-X)^(N-1)/[X^N+(1-X)^N]^2.
With |Y|<=sum_mu |Phi_mu|^2<=4 and the stated L2 jet norm,

    integral |Y|^N <=4^N U_Phi^2.

Therefore the complete nonconstant profile action is bounded
by [8epsilon 4^N/kappa0^(N-1)] integral U_Phi^2.

Its scalar force follows by differentiating the fixed
coefficient, with deltaY=2 partialPhi.partialeta.
The three contributions are bounded by

    [8epsilon 4^N/kappa0^(N-1/2)
     +64N epsilon 4^(N-1)/kappa0^(N-1)
     +16epsilon 4^N/kappa0^N] ||Phi||J3||eta||J3.

For the first term use ||Y^N||2<=4^N U_Phi; for the
gradient terms use the sup bounds on the remaining
powers of Y and
integral 2sum_mu |Phi_mu eta_mu|<=8U_Phi U_eta.
This controls the full function, not only its vacuum germ.

Exact inequalities 4^1024<10^617 and64N<10^5 give action
coefficient below10^-818552 and summed force coefficient
below10^-818547. Both are below the loose display1e-818500.
The report records exact exponent arithmetic and symbolic
anchor powers instead of printing million-digit denominators.

At fixed physical Minkowski g the canonical vector quadratic
operator and measure do not depend on Phi. Their fixed
renormalized determinant therefore has zero scalar variations,
including after regulator removal in this specified free
sector. Source/contact mean and scalar-force covariance
remain as computed; they are not removed with the determinant.

Along the actual family, all canonical scalar functions
f+Deltaf and r, the vector mass and its fixed prescription
stay fixed. Complete source mean force and noise variance
decay as1/kappa on this fixed finite-time class; this is a
conditional Gaussian decoupling result. Deltaf itself
PERSISTS in the decoupled scalar action. No interacting
light/graviton/auxiliary quantum limit, full scattering
contour or loop-truncation estimate is inferred.
