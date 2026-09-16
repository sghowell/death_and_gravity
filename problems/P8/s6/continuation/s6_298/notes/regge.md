# Whole-Q order error and a conditional transfer-contour budget

Write v=cosh eta>1. In conventional, not Olver-rescaled, normalization,

Q_nu(v)=integral_0^infinity
[v+sqrt(v^2-1)cosh y]^(-nu-1)dy.

This is [DLMF14.12.6](https://dlmf.nist.gov/14.12.E6) at order zero after
multiplying Olver's Q by Gamma(nu+1). Put
A=v+sqrt(v^2-1)cosh y=e^eta/u and d=e^-2eta. Its positive form is
Q_nu=e^[-(nu+1)eta] integral_0^1 u^nu/ sqrt[(1-u)(1-du)]du.

Consequently Q_(2+delta)/Q2 is the expectation of
exp[-delta(eta-log u)] for the normalized positive u^2 weight.
The mean of -log u is at most the mean of(1-u)/u. At d=0 the latter
is (4/15)/(16/15)=1/4. For0<d<1, the extra weight(1-du)^-1/2
increases with u while(1-u)/u decreases, so the covariance identity
makes this expectation no larger. Using1-e^-x<=x gives, for every delta>=0,

0<=1-Q_(2+delta)(v)/Q2(v)<=delta(acosh v+1/4).

This is a whole-function inequality. Its logarithm diverges at a pair
threshold; a uniformly small relative error there is not assumed.

Now CONDITIONALLY take one isolated even pole, continued two-channel
unitarity and alpha(T)=2+delta(T),0<=delta(T)<=T/L. Factor the nonkinematic
species-coupling ratio as zeta_a, with zeta_light=1 and
|zeta_heavy-1|<=epsilon. The pole-coupling kinematic ratio is
[(T-4a)/(T-4mu)]^(alpha/2). Relative to the exact spin2 cut its product
with the continued Bose wave is

zeta_a [(T-4a)/(T-4mu)]^(delta/2) Q_(2+delta)(v_a)/Q2(v_a).

The two latter factors lie in[0,1]. Its absolute defect is therefore at
most epsilon+delta L_a(T), where

L_a=acosh(v_a)+1/4+log[(T-4mu)/(T-4a)]/2.

Let ell(t) denote the selected one-loop correction to log[f(t)/f(0)].
This normalized definition retains the subtraction needed for the logarithmic
derivative; it is not a new finite renormalization condition. Assume its
slit transfer disk to U>4n contains only the stated cuts and has arc norm
at most B. No zeros or other singularities of the normalized residue may
be hidden in this premise. Cauchy's derivative formula and the exact
positive massive graph tail then imply

|ell'(0)-2f1'(0)| <=
B/U + g^2/(16pi^2 n U) + epsilon*g^2/(144pi^2 n^2)
+ g^2/(480pi^2 n^3 L) *
 {q^2[log(A/q)+3/4]+x^2[log(A/x)+log(U/x)/2+1]},

where q=U-4mu,x=U-4n,A=2U+4n.

To prove the last terms, both exact spin2 densities obey
2 Im f1_a/(pi T^2)<=g^2(T-4a)^2/(240pi^2 n^3 T^2);
for HH the q_a/q_phi factor cancels the q_phi in b^2.
Use acosh v<log(2v), bound its numerator by A, and for HH add
half log[U/(T-4n)]. Multiply delta<=T/L and use
(T-4a)^2/T<=T-4a. Integrating y log(A/y) gives
y^2 log(A/y)/2+y^2/4. The threshold endpoints vanish.
The displayed positive bound follows for the full intervals.

Dividing by kappa bounds only the selected finite contour-coefficient
mismatch: -ell'(0)/kappa versus the known -2f1'(0)/kappa.
It does not bound the other terms of the gravitational sum rule.

The original model has NOT supplied U,L,epsilon or B. The complex-J
continuation, residue ratios, absence of additional cuts, analytic trajectory
matching and the complex-energy high-energy remainder are not derived.
This conditional map cannot be used as a numerical P8 gravity verdict.

The motivation is the one-Regge-tower mechanism in
[Caron-Huot and Tokuda](https://arxiv.org/html/2406.07606v2).
Here the approximation error and omitted arc are explicitly retained;
the paper's light-mass leading-term approximation is not adopted as a bound.
