# A causal logarithmic inverse with a quantitative norm

This is an isolated dimensionless operator theorem. It supplies one tool
for an actual semiclassical fixed-point argument; it does not supply that
argument's nonlinear map, compatible data, interval or smoothness theorem.
The parameter below does not change the fixed A.8 stress prescription.

## 1. Operator and correct domain

For real `beta`, `f in C1[0,L]` and Euler's constant `gamma_E`, define

\[
 (D_\beta f)(t)=\frac1{8\pi^2}\left\{
  (\beta-\gamma_E)(f(t)-f(0))
  -\int_0^t\log(t-s) f'(s)\,ds\right\}.                 \tag{1}
\]

The integral is a usual locally integrable logarithmic convolution.
`D_beta f` is continuous and vanishes at zero. We do **not** assert that
it belongs to `C1`: for example

\[
 D_\beta[t]=\frac{t}{8\pi^2}
   (\beta-\gamma_E+1-\log t),                            \tag{2}
\]

whose derivative diverges at zero. In particular no `C1 -> C1`
endomorphism or automatic smoothness conclusion is used here.

Extending `f-f(0)` beyond `L` with exponential order, the Laplace transform
of (1) is

\[
 \widehat{D_\beta f}(s)
  =\frac{\beta+\log s}{8\pi^2}
       \left(\widehat f(s)-\frac{f(0)}s\right).           \tag{3}
\]

The logarithm is its principal branch. The identity follows from
`Laplace(log t)=-(gamma_E+log s)/s` and the derivative transform, including
the initial value. An extension beyond `L` cannot affect any retarded
convolution below `L`. The exact polynomial-data checks retain both the
harmonic numbers and Euler's constant, rather than suppressing a finite
term in (3).

## 2. Pole and cut, both required

Put `p=exp(-beta)>0`. The causal inverse kernel has the positive form

\[
 K_\beta(t)=8\pi^2\left[p e^{pt}
  +\int_0^\infty\frac{e^{-rt}}{(\log r+\beta)^2+\pi^2}\,dr\right],
 \qquad t>0.                                             \tag{4}
\]

One way to prove (4) without treating a conditionally convergent Fourier
integral as absolutely convergent is the keyhole identity

\[
 \frac1{\beta+\log s}
   =\frac p{s-p}
       +\int_0^\infty
          \frac{dr}{(s+r)((\log r+\beta)^2+\pi^2)},
 \qquad s>p.                                             \tag{5}
\]

Apply the residue theorem to the principal-log function divided by
`z-s`, using circles of radii tending to zero and infinity and the two
banks of the negative real cut. The outer integral is `O(1/log R)`;
the inner integral is `O(epsilon/abs(log epsilon))`. The pole at `p`
has residue `p`, and the cut values `log r +/- i*pi` give density
`1/((log r+beta)^2+pi^2)` with the positive sign in (5).
Small circles around `s` and `p` give the displayed two residues.
All integrals in (5) converge. Taking Laplace transforms of the two
nonnegative terms in (4), Tonelli's theorem gives (5) directly.

Thus, with

\[
 (J_\beta h)(t)=\int_0^t K_\beta(t-s)h(s)\,ds,             \tag{6}
\]

we have `J_beta D_beta f=f-f(0)` for `C1` data, by (3), (5) and Laplace
uniqueness. Causality makes this a local-interval statement. For general
continuous `h`, (6) is the inverse in the causal-distribution sense;
it need not be a `C1` solution of (1).

In particular omitting the positive pole would leave a Laplace defect
`p/(s-p)`. It cannot be removed while retaining the causal inverse (3).
This also gives the lower norm bound
`8*pi^2*(exp(p*L)-1)`: a small-interval inverse estimate is not a
uniform-in-duration stability estimate.

## 3. Exact norm and finite bound

Positivity gives the exact operator norm on `C[0,L]`:

\[
 \|J_\beta\|_{C^0\to C^0}
  =C_\beta(L)=8\pi^2\left[e^{pL}-1
    +\int_0^\infty\frac{1-e^{-rL}}
       {r((\log r+\beta)^2+\pi^2)}\,dr\right].            \tag{7}
\]

The upper bound follows by convolution with `abs(h)` and is attained
by `h=1` at `t=L`. On the closed subspace `h(0)=0`, approximate that
function from below away from zero to obtain the same supremum norm.

For `beta-(log L)/2>0`, split the cut integral at `R=L^(-1/2)`.
Below `R`, use `1-exp(-rL)<=rL` and the denominator at least `pi^2`.
Above `R`, use `1-exp(-rL)<=1` and the denominator at least
`(log r+beta)^2`; its antiderivative is `-1/(log r+beta)`. Hence

\[
 C_\beta(L)\le 8\pi^2\left[
 e^{pL}-1+\frac{\sqrt L}{\pi^2}
       +\frac1{\beta-\tfrac12\log L}\right].             \tag{8}
\]

This proves local integrability of the kernel, boundedness and continuity
of (6), and `C_beta(L)->0` as `L->0`. Integrability on any larger compact
interval then follows directly from (4) away from zero.

There is no automatic derivative gain of one order. For example the cut
part of (4), integrated over `t/2<=s<=t` and `1/t<=r<=2/t`, bounds
`J_beta[1](t)` below by a positive constant divided by
`(abs(log t)+abs(beta)+log 2)^2+pi^2` for small positive `t`.
Consequently `J_beta[1](t)/t` is unbounded and this output is not `C1`
at zero. This is a domain control, not an obstruction to a separate
smooth-data theorem with compatible initial conditions.

## 4. Fully rational dyadic implementation

For integer `B>=0`, `beta>=-B`, and integer `n>2B`, set
`L=2^(-2n)` and `E=3^B L`. If `E<1`, (8) is bounded by the rational

\[
 C_{B,n}^{\rm rat}
    =\frac{80E}{1-E}+\frac8{2^n}+\frac{80}{n/2-B}.        \tag{9}
\]

Here `exp(-beta)<=exp(B)<=3^B`, `exp(E)-1<=E/(1-E)`,
`log 2>1/2` and `pi^2<10`. The first two inequalities follow from
positive exponential series and `j!>=1`; `e<3` follows by comparison
of its tail to a geometric series. The logarithm inequality follows
by integrating `1/x>=1/2` on `[1,2]`. For the last inequality the
positive integral

\[
 \int_0^1\frac{x^4(1-x)^4}{1+x^2}\,dx=\frac{22}7-\pi>0,
 \qquad 10-(22/7)^2=6/49>0,
\]

is checked by exact polynomial division in addition to symbolic
integration. No rounded floating-point constants enter (9).

The deliberately conservative example `B=1,n=512` has
`L=2^(-1024)` and `C_rat<1/3`. It includes `beta=gamma_E`, since
`0<gamma_E<1`. This very short **dimensionless** interval is only an
operator calibration; it is neither an actual proper-time duration
nor a calibrated physical semiclassical evolution interval.

## 5. Conditional a-posteriori theorem

Let `X=C_0[0,L]`, with a closed radius-`r>0` ball about `xbar in X`.
Suppose a concrete map `F` has been proved to map this ball to `C[0,L]`
with uniform-norm Lipschitz constant `Lambda`, and set `C[x]=J_beta F[x]`.
If independently established inputs satisfy

\[
 q=C_{\rm upper}\Lambda<1,\qquad
 \eta=\|C[\bar x]-\bar x\|_\infty,\qquad
 \eta+qr\le r,                                           \tag{10}
\]

then the map preserves the ball and Banach's theorem gives a unique
fixed point there, with

\[
 \|x-\bar x\|_\infty\le\frac\eta{1-q}.                  \tag{11}
\]

The code checks these rational inequalities and rejects `q>=1` or a
failed self-map condition. It does not verify a user's unproved
Lipschitz or residual input. The illustrative values
`C_upper=1/3,Lambda=1,eta=1/10,r=1` give distance at most `3/20`.
They are an abstract example, not inputs derived for the SEE map.

For the physical problem one must still construct the full actual-state
map with its anomaly and constraint, propagate compatible metric/state
data, prove the required regularity, and enclose the desired interval.
The companion actual-mode lemma controls one nonlinear part of that map;
it does not bound the remaining, potentially stiff gravitational terms.
