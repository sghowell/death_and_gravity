# Smooth support and a quantified finite-profile budget

Use exactly the S6.55 C-infinity step chi(s), zero for s<=1
and one for s>=2, with the flat exponential transition

    chi(s)=b(s-1)/[b(s-1)+b(2-s)], b(t)=exp(-1/t), t>0.

Define w(x)=1-chi(64(x+1)^2). The polynomial argument avoids
an absolute-value corner. The window is one for |x+1|<=1/8
and zero for |x+1|>=sqrt(2)/8. Thus the original tube
|x+1|<=1/10 lies inside an open flat region. An open
neighborhood of the Minkowski clock norm x=0 is untouched,
including every local field jet there. This preservation of
vacuum coefficients is not a vacuum-existence, common-parent or
UV-matching proof.

## Finite, continuous window derivatives

S6.55 supplies b^(j)(t)=exp(-y)Q_j(y), y=1/t. For each monomial
y^r, exp(-y)y^r<=r! follows from the positive exponential series.
Thus B_j=sum |q_{j,r}|r! bounds the j-th bump derivative.
In the transition one of s-1,2-s is at least 1/2, so the
denominator is at least exp(-2)>1/8. The last strict inequality
is checked by the exponential series through order six, bounding
its remaining tail by a geometric series of ratio 1/4.

For reciprocal-denominator derivative envelopes set R_0=8 and

    R_j=16 sum_{k=1}^j binomial(j,k) B_k R_(j-k).

Differentiating numerator times reciprocal gives step envelopes
C_j=sum_{k=0}^j binomial(j,k)B_k R_(j-k), j>=1.
Use C_0=1, since 0<=chi<=1. On the transition |x+1|<1/4,
so the inner polynomial has first derivative bounded by 32,
second derivative 128 and no higher derivatives. The exact
quadratic composition rule gives

    W_j=sum_{k=0}^{floor(j/2)}
          j! C_(j-k) 32^(j-2k)64^k/[(j-2k)!k!], W_0=1.

This rule is independently checked through order five. The
resulting global derivative envelopes for w are

    1, 8704, 22874112, 92896886784,
    509829585567744, 3527551905599324160.

They are conservative bounds over the continuum; derivatives
at the transition endpoints vanish to every order.

## Mixed profile bounds

Let eta_j bound both |rho_sigma^(j)| and |p_sigma^(j)| from
S6.59. Inside the clock tube the dimensionless profile obeys

    |partial_u^j P|<=11eta_j/10,
    |partial_u^j P_x|<=eta_j, P_xx=0.

On the entire support |x+1|<1/4, its affine inner factor is
bounded by 5eta_j/4 and its first x derivative by eta_j.
The product rule therefore gives, for j+k<=5,

    |partial_u^j partial_x^k P|
       <=[(5/4)W_k+k W_(k-1)]eta_j,

where the second term is absent for k=0. Outside the support
the profile and all derivatives vanish. At L=10^24,m=1000
all 21 bounds are below 10^-18, including the smooth
off-tube transition. The largest is below 1.456*10^-19.
The exact rational bounds, not these displays, are certified.
These are coordinate coefficient-jet bounds, not canonical
interaction or Wilsonian-cutoff estimates.

For the actual point-chart lapse-square coefficient put y=1/h
in [0,1]. Its energy coefficient (3y-1)/2 has absolute value
at most one. Its pressure coefficient is -(21y^2-18y+4)/8.
The identities

    21y^2-18y+4=21(y-3/7)^2+1/7,
    7-(21y^2-18y+4)=(1-y)(21y+3)

give absolute value at most 7/8. Thus |Delta J|<=15eta0/8.
The nv and v^2 coefficients are bounded by 15eta0/2 and
9eta0/2. Since the predecessor J_e>3/8 on the window,
its literal lower-scalar lapse-square coefficient stays above
1/3 in the scale example. This bound adds no assertion about
the full quantum constraint system or its causal cone.
