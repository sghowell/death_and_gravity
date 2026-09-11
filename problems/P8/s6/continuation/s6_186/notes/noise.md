# Actual centered stress covariance and its normalization

The unchanged all-order CD state is the pure quasifree
state determined by its prepared positive mode data.
For its field expansion in creation and annihilation
modes, the centered local quadratic stress applied to
the state has only the two-particle creation part.
With the Fourier convention in notes/oscillatory.md,

    Var T[f]=2 sum_(r,s=1)^3 integral d^3k d^3l/(2pi)^6
      |integral dt a^3 sum_ab fhat_ab(t,k+l)
        conjugate(u_r(t,k))^T M_ab conjugate(u_s(t,l))|^2.

The factor2 is the pair-exchange contraction. All nine
polarization pairs and all tensor cross terms are
retained. The conjugate mode amplitudes are essential;
the phase sign does not alter the oscillatory bound.

Split each pair into reference and exact remainder
and use |x+y|^2<=2|x|^2+2|y|^2. The total prefactor
is2*9*2=36. Complete internal integrals and spatial
Plancherel give

    Var T[f] <=36*1e50/m sum_(j=0)^3||partial_t^j f||L2^2
       +288*1e30/m^7 (||f||L2^2+||grad f||L2^2).

At m1000 this is below1e50 N[f]^2, where

    N[f]^2=sum_(j=0)^3||partial_t^j f||L2(dt dx;F)^2
              +||grad_spatial f||L2(dt dx;F)^2.

Thus STD(T[f])<1e25 N[f]. Both normalizations are
explicit: STD(T[f]/kappa)<1e-775 N[f] and the
Einstein-normalized metric test coupling has
STD<1e-375 N[h]. The latter uses
delta g_ab=2h_ab/sqrt(kappa) together with the
variational factor1/2. It is not yet a fully
constraint-reduced mixed-mode canonical norm.

The two-particle kernel is square integrable by these
bounds. Regulated momentum integrals therefore converge
to this centered smeared Wick observable. Tests are
smeared before the limit. No pointwise stress variance,
finite momentum cutoff or product of independently
smoothed fields replaces the local stress.

The same fixed covariant stress prescription remains.
At a fixed metric, the local subtraction/finite terms
and the fixed scalar retuning are c-numbers; they
cancel in centering. Using actual-state Wick ordering
to compute the centered operator does not redefine
its prescribed one-point mean.

For context, [Hu and Verdaguer, arXiv:0802.0658v1,
section3.2](https://arxiv.org/pdf/0802.0658) distinguish
centered symmetrized stress covariance from response
and metric dynamics and emphasize its distributional
meaning. Their scalar kernel is not imported here.
The actual Proca tensor and its numerical constants
are derived from the constrained modes above.
