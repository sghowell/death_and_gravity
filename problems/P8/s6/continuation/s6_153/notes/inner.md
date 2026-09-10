# Fixed physical inner subtraction in d dimensions

At the scalar physical mass one, put

    Delta(s,x)=xM+1-x-x(1-x)s,
    D(x)=Delta(1,x)=xM+(1-x)^2,
    b(x)=x(1-x)/D(x), ell=log(mu^2).

The full mixed scalar bubble is (g/Q) exp(gamma e+ell e)
Gamma(e) integral Delta^(-e) dx, Q=16pi^2.
Differentiate in s BEFORE removing the regulator:

    alpha_D=Pi_D'(1)
      =(g/Q) exp(gamma e+ell e) Gamma(1+e)
        integral b D^(-e) dx.

This is the asymptotic multiplier in the fixed physical
on-shell-subtracted kernel Pi_R(s)/(1-s).
The gamma factor exp(gamma e)Gamma(1+e) has zero linear
coefficient. Thus

    alpha_0=(g/Q) integral b dx,
    alpha_1=(g/Q) integral b[ell-logD] dx.

All denominators are strictly positive on the compact parameter
interval, so these derivatives are uniform near e=0.
They use precisely the parent S6.112 parameter weight.

For large Euclidean t=-s, the on-shell quotient is
alpha_D+[Pi_D(-t)-Pi_D(1)]/(t+1).
The second term decays. Uniformly for epsilon in a sufficiently
small compact neighborhood of zero, the massive Feynman
parameter integral grows at most as a fixed power t^delta
with delta<1; a possible logarithm at zero epsilon is covered
by the same bound after the constant UV pole is cancelled.
The quotient remainder is O(t^(-1+delta)) plus a constant/t
term. Its insertion in the outer two-propagator loop is
uniformly integrable when the epsilon neighborhood is small
enough, for example |Re e|<1/8. Thus only alpha_D multiplies
the remaining logarithmic outer reference.

The inner tadpoles and all affine mass/kinetic counterterms
are annihilated by Pi(s)-Pi(1)-(s-1)Pi'(1) before estimating.
No physical mass, residue or heavy-source choice changes here.
