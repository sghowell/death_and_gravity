# Complete global exponential type at the exact matter radius

On either regular segment, set V=[1,0;-sign*w/d,1] and
omega=diag(c/a,1/a), c^2=J/(J+delta_J). The explicit mode matrix
S=[V,V;i*V*omega,-i*V*omega] diagonalizes the principal generator of
Y=(k*x,x'). Its determinant is -4*c/a^2, with no inverse frequency-gap
factor. The global lower c^2>=1199/1215>0 and c<1 hold at every finite
time. Therefore S,S^-1,S' have finite norms on each compact regular
segment, including where the clock and matter speeds approach one
another asymptotically. No assertion at infinite time is needed.

The COMPLETE scaled velocity generator is

    Y'=[k*J0+R(k,u)]Y,
    J0=[0,I;-K^-1*G/a^2,0].

R has bounded complex-frequency norm for |k|>=64: its lower blocks are
-(Q-q*K^-1*G)/k and -L. The rational-degree checks, explicit pole
separation, and compact smooth coefficient domains prove this bound
uniformly, not merely along positive real k. Compactifying 1/k adds
the finite limits at infinity; no pole occurs in that closed domain.

In mode coordinates the generator is i*k*diag(c/a,1/a,-c/a,-1/a)
plus a bounded remainder. Its Hermitian principal part is
-Im(k)*diag(c/a,1/a,-c/a,-1/a). Euclidean logarithmic-norm Gronwall
therefore bounds each COMPLETE segment transfer by a finite compact-strip
constant times exp(|Im(k)|*integral du/a). The real part of k causes no
exponential growth. Endpoint canonical maps add at most three powers
of |k| per segment. At most three segments give

    ||U_rho(t,s,k)|| <= C_R*(1+|k|)^9*exp(S_m(t,s)*|Im(k)|), |k|>=64.

The segment radii ADD to the exact S_m=F(t)-F(s), because the evolution
is composed at the same actual times. For |k|<=64 the original polynomial
generator has a finite norm on the whole compact strip and gives a
uniform bound. This low-root region includes arbitrarily large complex
null vectors and every finite-q velocity-chart pole.

For k^2=zeta dot zeta, the pinned exact Gram identity gives
|k|<=|zeta| and |Im(k)|<=|Im(zeta)|, with either pointwise root.
Consequently the entire original scalar multiplier has global bound

    |G_hat(t,s,zeta)| <= C_R'*(1+|zeta|)^9*exp(S_m(t,s)*|Im(zeta)|).

Apply S6.97's pinned entire-type-to-distributional-support argument:
the inverse spatial Fourier distribution is supported in the radius-S_m
ball. Its uniform polynomial real-momentum bound and continuous time
dependence define a joint retarded spacetime distribution. Ordered-time
support, zero initial chi-to-pi_chi transfer, and time antisymmetry give
the stated matter-causal response. Every finite pair lies in some finite
strip, so this proves the global-clock statement without a uniform bound
as |u| tends to infinity.

F'=1/a is checked exactly and F(+infinity)-F(-infinity)=pi/2. This finite
conformal-time length is compatible with the already established complete
physical classical bounce; it is not a finite proper-time endpoint.
