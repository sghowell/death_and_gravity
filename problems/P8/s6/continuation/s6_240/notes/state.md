# Full exact heavy state, physical phase and reference readout

## Reference equation and physical normalization

On a=(1+t²)², Hubble=4t/(1+t²), the complete minimally coupled mode equation is

S_p''+3Hubble S_p'+omega² S_p=0, omega²=n+p²/a².

For chi=a^(3/2)S this becomes chi''+(omega²-U)chi=0, U=3Hubble'/2+9Hubble²/4. Globally Hubble²<=4 and Hubble'<=4, so U<=15<n. The metric is smooth, nondegenerate and globally hyperbolic at every finite real t, and n>0. All p>=0, including zero momentum, are retained.

Physical positive frequency is exp(-i omega t). The normalized physical Wronskian is S conjugate(S')-conjugate(S) S'=i/a³, with W2(t,t')=S(t)conjugate(S(t')). For Cauchy vector (S,a³S'), its positive rank-one Gram matrix has antisymmetric part i times the canonical symplectic matrix. These signs are independently checked and not taken from the opposite mode convention in a literature formula.

The physical mode readouts, derived from the full lapse and scale-factor variations, are

rho_p=(|S'|²+(n+p²/a²)|S|²)/2,

P_p=(|S'|²-(n+p²/(3a²))|S|²)/2.

Their exact equation implies rho_p'+3Hubble(rho_p+P_p)=0. They are integrated with p²dp/(2pi²). The zero reference coherent mean does not remove fluctuations or the metric response of the determinant.

## Exact compactly smeared minimizer

Fix x=8(t+5/8), phi=exp[-1/(1-x²)] on |x|<1 and0 outside, w=phi²/Z, Z=integral phi²dt. This gives a real C-infinity sampling amplitude f=phi/sqrt(Z), normalized w=f² and support[-3/4,-1/2]. On |x|<=1/2, exp[-2/(1-x²)]>exp(-3)>1/27, hence Z>1/216.

For ANY exact normalized comparison basis S define

c1=integral w(|S'|²+omega²|S|²)/2,

c2=integral w(S'^2+omega²S²)/2.

Writing S=x+iy,S'=v+iw gives the exact cone identity
e²-|q|²=omega²(xw-yv)²=omega²/(4a^6)>0.
The future Lorentz cone is convex; integrating the strict inequality against a nonzero positive weight gives c1>|c2|. More directly c1>=integral w omega/(2a³)>=Omega_star/128 on[-1,1], where Omega_star=sqrt(n+p²/16), since a<=4.

Set delta=sqrt(c1²-|c2|²),
alpha=sqrt((c1+delta)/(2delta)),
beta=-c2/sqrt[2delta(c1+delta)],
and T=alpha S+beta conjugate(S).
This is regular at c2=0. The full transformation identities give alpha²-|beta|²=1, c2[T]=0 and c1[T]=delta. Any further normalized Bogoliubov transform of T has energy delta(1+2|b|²), so this is the global minimizer, unique up to an overall phase. Replacing the initial exact basis therefore leaves the covariance unchanged.

The hypotheses of Olbermann's existence/minimization and Hadamard theorems apply: a smooth Robertson-Walker geometry, a positive massive minimally coupled real scalar, and a nonzero real compactly supported smooth sampling amplitude. This application to the displayed metric and sampling is our inference after checking those hypotheses. The theorem supplies Hadamard regularity of the full exact selected state. It supplies none of this package's numerical constants. See [the primary source](https://arxiv.org/abs/0704.2986).

## Fixed family, source and scope

Select this reference state ONCE. On smooth compact perturbations of the metric/source strictly to the future of t0=-1/2, sharing its Cauchy neighborhood and global hyperbolicity, keep the exact reference Cauchy covariance. Propagate the homogeneous KG equation and add the unique retarded smooth coherent solution for the complete heavy source. Causal propagation preserves CCR and positivity; the propagated Hadamard singularity and the added smooth one-point function have the usual fixed-state meaning. No new minimization is performed on the live history.

The S238 classical source g Phi²(1-X)^8 exp(-A X²)/2 and the S239 finite onepoint switch vanish at the clock with the required jets. Thus the reference heavy mean and direct source-induced clock force are zero. Classical coherent mean begins at history order eight; eliminating it gives the first source-induced Euler force at order fifteen. The free heavy action and its connected metric response remain.

WKB iterates in the bound proof are only a comparison basis for exact modes. A finite WKB state is never identified with this selected Hadamard state. The definition is global on the reference; the quantitative bound is restricted to[-1,1]. The full interacting light/metric quantum state remains open.
