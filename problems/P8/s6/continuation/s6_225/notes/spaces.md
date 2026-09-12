# Uniform Sobolev bounds, causal graph and physical density

Fix a conformal slab of lengthT<=1 and r real. Use the two scalar amplitudes(w,c) and their dual Euclidean source norm at each Fourier momentum. Spatial Fourier transformation is unitary after the chosen Sobolev weight; q=P^2 and the Pi embedding is bounded uniformly in direction. The exact spatial metric Gram matrix is[[12,-4],[-4,4]], between2I and14I. Thus converting amplitude and dual-source estimates to other physical Frobenius conventions requires the retained embedding constants; none are silently identified with1.

The scalar dominating time kernel for Ecurv is

e(tau)=375 tau^3 exp(108tau)/16, 0<=tau<=T.

Minkowski in time and Plancherel in space give

||Ecurv||_(C_tH^r to C_tH^r)
 <=integral_0^T e(tau)dtau
 <=375T^4 exp(108T)/64 <6T^4 exp(108T).

Young's inequality for the same dominating kernel gives the identical L2_tH^r bound. This is uniform in ALL comoving momenta with no spatial derivative loss, on these specified amplitude spaces. Continuous dependence in time follows first on compact momentum sets and then by domination. Differentiating the outerY factor also gives a finite first-time-derivative bound becauseY(t,t)=0 and partial_tY is uniformly bounded.

## Complete initial-boundary forward graph

At each fixedq, use causal distributions with no forgotten initial atoms. Y andZ are the unique smooth-coefficient causal ODE inverses. The scalar forward factors are the explicit S224 causal distributions

Fi(D^2+q)=-Ai(0)delta-(D^2+q)Gi,q,

where Gi,q is uniformly bounded and continuous:9pi/(128m) in shear and27pi/(64m) in trace. The distributional derivatives retain the initial boundary and avoid a divergent split of local/cut pieces. Their inverse isD J_i,q with the proven bounded primitive.

All these factors act on time distributions with values in spatial Sobolev spaces, allowing a finite reduction of the spatial index when derivatives are applied. Smooth coefficient multiplication and causal integration preserve this union of finite-order Sobolev-valued distributions. Y/Z's higher finite time jets have at most polynomialq growth by their ODEs. The scalar forward factors have at most one explicitq factor in the displayed representation, while the coordinate factors each have at most oneq factor. In particular Aref sends continuous H^r amplitudes to finite-order time distributions with values inH^(r-6), a deliberately crude finite-order bound.

The graph domain consists of continuous zero-past amplitudesu with u(0)=0 such that the complete causal distribution Aref u is an ordinary sourcef inC_tH^r, including the initial boundary. The ordinary Ecurv integral maps every such source to a continuous zero-past response. The distributional factor identities give Aref Ecurv f=f. Conversely Ecurv Aref u=u proves uniqueness on the graph. No density or invariance of the FULL S222 quantum-force graph is inferred. One may first prove these identities for smooth compact-momentum inputs and pass in the finite-order distribution topologies; the no-loss norm applies to the resulting ordinary reference inverse.

ForL2 sources the same inverse and identities are understood in the corresponding distribution graph, with the continuous representative supplied by its integral kernel. This is not a claim that Aref is a bounded forward self-map onC_tH^r.

## Physical conformal density is not commuted

The reference force normalization in conformal measure is defined as

Qbar_ref=M_[1/(64pi^2 kappa a^4)] Aref.

Its inverse is Ecurv M_[64pi^2 kappa a^4], on the RIGHT. Indeed Qbar_ref Ecurv M=I and Ecurv M Qbar_ref=I by adjacent multiplication, not by moving M through Ecurv. Sincea^4<6 on the unchanged slab, the norm is below

2250pi^2 kappa T^4 exp(108T).

The exact matrix controls detect the reversed multiplication order. Kappa remains10^800, so this is not a small-feedback estimate. Full force/chart contacts and the actual curved nonlocal factors must still be matched to S222 before using this reference inside that system.
