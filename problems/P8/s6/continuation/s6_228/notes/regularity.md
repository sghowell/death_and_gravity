# Smooth preparation, retained boundaries and the original equation

The bounded Neumann construction first gives a continuous solutionh of

h+K2Vh=K2 f0, f0=I4g.

Causality preserves the common nonempty zero initial neighborhood ofg: every Neumann term vanishes there. This fact is used as an original preparation condition, not a reset of arbitrary higher-derivative data.

Write V_j=(Dt+Ds)^j V. The proved boundsCj(1+|log(t-s)|) make eachV_j a bounded operator on continuous functions over the compact interval. The simultaneous derivative is essential. Differentiating the singular factorlog(t-s) in only one variable would create a nonintegrable1/(t-s) kernel, and no such absolute estimate is used.

## Time-shift proof of differentiability

LetS_epsilon h(t)=h(t+epsilon), initially on the shortened output interval. For epsilon smaller than the original initial zero neighborhood, changing the source variable gives

S_epsilon(Vh)=V_epsilon S_epsilon h,
V_epsilon(t,s)=V(t+epsilon,s+epsilon).

There is no lost initial term: the small extra source interval after the change of variable lies entirely in the retained zero germ. The convolutionK2 commutes with the same prepared shift. Subtract the shifted and unshifted equations, keeping the originalV on the left:

(I+K2V) Delta_epsilon h
 =K2 Delta_epsilon f0
  -K2[(V_epsilon-V)/epsilon] S_epsilon h.

Away from the diagonal, the difference kernel converges uniformly on compact sets toV_1. Near the diagonal it is dominated by the integrable weak-log majorant supplied by the simultaneous derivative bound. Splitting the triangle into these two regions gives convergence in the kernel operator norm on continuous functions. S_epsilon h tends uniformly toh and the inverse ofI+K2V has the same uniform finite bound on every shortened interval. Hence Delta_epsilon h converges uniformly to a continuous derivative. The one-sided limit supplies the endpoint derivative without a final-time condition.

Repeat this argument. The derivative rule is

D^n(Vh)=sum_(j=0)^n binom(n,j) V_j h^(n-j),

and prepared convolution byK2 commutes with ordinary time differentiation. Therefore

(I+K2V) h^(n)
 =K2 f0^(n)
  -K2 sum_(j=1)^n binom(n,j) V_j h^(n-j).

Induction gives every fixed finite time derivative, with finite bounds depending on the correspondingCj and source derivatives. No uniform-in-order analytic norm or numericalCj is asserted. Thus smooth preparedg gives smooth preparedh.

An independent weak-log exampleV(t,s)=(1+t+s)log(t-s) checks the simultaneous derivative rule. A second test with nonzero initial source value retains and detects the otherwise missing termV(t,0)h(0). The latter is a boundary diagnostic, not an enlargement of the admitted prepared problem.

## Both identities for the actual unintegrated operator

For smooth preparedh, the exact physical normal-form identity proved from the current is

I4 T_total h=(F2+V_total)h.

The constructed smoothh obeys(F2+V_total)h=I4g. ApplyingD^4 to the full causal distributions givesT_total h=g. Zero-past I4 andD^4 retain all initial distributions; they are not a projection onto a subset of higher-order modes.

Conversely, any smooth preparedh satisfying the original equation gives the normal-form equation afterI4. The original causal identityK2F2=I and the bounded inverse ofI+K2V_total then identifyh with the constructed solution. This proves existence, uniqueness and both inverse identities on the stated smooth prepared homogeneous class. No arbitrary final condition, order reduction, discarded pole, or separate choice of extra-branch data enters.

The five tracefree components square-sum in any fixed orthonormal Frobenius basis and share the same kernel. This remains a literal zero-transfer time problem per comoving volume. It does not establish a uniform spatial Fourier multiplier, the maximal distributional domain of the full coupled system, or S222 graph invariance.
