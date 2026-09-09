# Actual coupled scalar energy and the two-chart cover

Use physical u=t/tau and z=1/q, where q is squared
physical momentum in the tau unit. Derivatives hold
comoving momentum fixed:

D_u=partial_u+2H z partial_z.

In the two scalar charts, the exact phase Hamiltonian
is the S6.76 scalar block after its two canonical
matter-boundary shifts. Its Legendre transform gives

L=a^3[dot(Q)^T alpha dot(Q)/2
      +dot(Q)^T beta Q-Q^T gamma_positive Q/2].

All coefficients are derived from the actual Hamiltonian
with J_e=J+4*10^-6/h^2 and w=-ell*Lambda. The gamma
canonical swap keeps the -H*b*P_b generator. No time
jet is frozen away.

Set S=(beta+beta^T)/2, A=(beta-beta^T)/2 and

V=gamma_positive+S'+3H*S.

After the symmetric time boundary, the exact equations are

alpha Q''+(alpha'+3H alpha+2A)Q'
 +(V+A'+3H A)Q=0.

The antisymmetric A is not dropped or replaced by a
diagonal oscillator approximation.

## High-frequency decomposition

Exact rational algebra gives V=qG+R and
beta=q beta0+beta_remainder. The leading beta0 is
symmetric; R and A have no high-frequency pole.
Their coefficients and required time derivatives are
rational functions of the background jets and z.

On the actual background,

G=[[2J/d^2+c^2,c],[c,1]],

with (d,c)=(Theta,-ell*Lambda/Theta) in the unitary
chart and (Lambda,ell) in the gamma chart.
The actual scalar background identities establish this
formula at all u, including the margin's distinction
between J in G and J_e in alpha.

The unitary kinetic matrix has the same square
decomposition with J replaced by J_e.
The gamma kinetic matrix equals its positive z=0
limit plus the exact positive rank-one term

z/D * vv^T,
D=Lambda^2-(J_e+ell^2 Lambda^2/2)z,
v=(sqrt(2)*(J_e+ell^2 Lambda^2/2)/Lambda,
   ell*Lambda/sqrt(2)).

## Continuous chart and denominator bounds

On I, |H|<=2, |ell|<=1/10, |Theta|<=3,
|Lambda|<=1/2 and 1/40<J<=J_e<50.
The gamma chart uses |u|<=1/4, where

|Lambda|>=1231/4913>1/4.

The unitary chart uses |u|>=9/40, where

|Theta|>=864/3125>1/4.

For z<=10^-4, D and Lambda^2-J_e*z exceed 1/32
on the gamma chart. These denominators are not
assigned such a bound on the unitary chart.
The remaining denominators are explicitly factored
into admissible chart factors and positive constants;
an unrecognized or zero-frequency denominator fails.

Each actual background jet through order two has an
exact positive-even denominator polynomial in u.
Absolute numerator coefficient sums on |u|<=1/2
give explicit rational bounds. Those are propagated
through the generic rational coefficient numerators
and certified denominator factors. No time samples
replace the continuous bounds.

Completed squares give

alpha>=I/1000, G>=I/1000,
||alpha||,||G||<10000.

Indeed the triangular square change has inverse
Frobenius norm squared at most three, while
2J/d^2>1/180. The eigenvalue bound 1/540 is
stronger than the retained 1/1000.

## Exact energy identity

For the real system define

E=a^3[Q'^T alpha Q'+Q^T V Q]/2.

The same identity holds for complex columns using
real parts of Hermitian products. Literal substitution
of the full coupled equations gives

E'=a^3[-Q'^T(alpha'+3H alpha)Q'
       +Q^T(V'+3H V)Q
       -2Q'^T(A'+3H A)Q]/2.

The instantaneous antisymmetric velocity term cancels
only in Q'^T A Q'; its derivative mixing with Q remains.

Since q'=-2Hq,

V'+3HV=q(G'+HG)+(R'+3HR).

Let the certified operator bounds for R,
alpha'+3H alpha, G'+HG, R'+3HR and A'+3H A be
B_R,B_a,B_G,B_r,B_m. Define

q_min=max(10^8,2000 B_R).

Then V>=qI/2000, and for q>=q_min the energy estimate is

|E'|/E <=1000 B_a+2000(B_G+B_r/q_min)+B_m/5.

The last factor uses sqrt(q)>=10^4 and
sqrt(1000*2000)<2000. Exact rational calculations
put this expression below 10^9 in both charts.
The upper is conservative, not an actual growth rate
or instability prediction.

Choose a half-window 1/(4*10^9). Energy then grows
by less than exp(1/2)<2. The latter inequality follows
directly by comparison of the exponential series with
the geometric series, without floating evaluation.
The scale-factor drift is also below a factor two.

The chart is selected at the window center:
gamma for |u_center|<=19/80, unitary otherwise.
A half-width below 1/160 keeps the whole window in its
respective domain. Only windows contained in I are used.
