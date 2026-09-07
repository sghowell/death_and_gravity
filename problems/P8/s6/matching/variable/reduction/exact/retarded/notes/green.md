# Finite-delta retarded history in the exact own-f tensor equation

This S6.34 result concerns the exact **linear** own-f tensor equation on
the S6.33 isotropic branch, with a prescribed off-shell physical-metric
anisotropy. It is not a freely evolving full-parent background, a response
to the free-chi matter source, or a fixed physical-frequency EFT theorem.
The equation and zero-data prescription are used directly; a retarded
kernel is not substituted into a single-copy action and declared to be a
variational effective action.

## 1. Exact normalization and admitted finite parameters

The norm-two tensor action, with primes denoting physical-g proper time,
is

\[
 L_T=M^2q'^2/4+K_fQ'^2/4-\beta_1b(Q-q)^2/2,
 \qquad K_f=M^2b^3/N.
\]

Here `q` is prescribed; `Q` is the hidden tensor amplitude. The own-f
equation is `(K_f Q')'+2beta1*b*(Q-q)=0`. Retain the full derivative of
`K_f`, without freezing it. Define

\[
 x=\frac{T}{\tau\sqrt\delta},\quad
 k=\frac{b^3}{N},\quad
 s=\frac{2Q_{\rm profile}b}{D/\delta},\quad
 0<\delta\le\frac1{625},\qquad -\frac14\le x\le0.                 \tag{1}
\]

The numerator `Q_profile` is distinct from the hidden tensor `Q`.
Since `beta1=(M^2/tau^2)*Q_profile/D`, multiplying the physical equation
by `tau^2 delta/M^2` gives exactly

\[
 (kQ_x)_x+s(Q-q)=0.                                                \tag{2}
\]

The independent Fraction proof in [bounds.md](bounds.md) proves on the
whole finite-parameter domain

\[
 \frac{19}{5}<k<\frac{21}{5},\qquad 42<s<65,
 \qquad 1\le D/\delta\le\frac32.                                 \tag{3}
\]

Indeed `v=delta*x^2<=1/10000`, so the exact stationary box applies. These
are finite-delta bounds, not an approximation by an inner limiting
operator. No numerical evolution or additional bound on `k_x` is used
below. The exact coefficient branch is smooth for each actual positive
delta; `delta=0` is never inserted as an action.

## 2. Retarded normalization and the first-zero argument

For a source time `a` in `[-1/4,0]`, let `R(x,a)` be the homogeneous
Sturm Green solution to the right of `a`:

\[
 (kR_x)_x+sR=0,\qquad R(a,a)=0,\qquad k(a)R_x(a,a)=1.              \tag{4}
\]

The jump is in the **flux**, not an unweighted derivative. Equivalently,
the regular first-order system is `R_x=V/k`, `V_x=-sR`, with initial data
`(R,V)=(0,1)`. Positive continuous `k` and continuous `s` suffice for this
system; the actual coefficients have more regularity. Write

\[
 k_-=19/5,\quad k_+=21/5,\quad s_+=65,\quad h=x-a\le1/4.
\]

Suppose the flux had a first zero in this interval. Before that zero,
`V>0`, so `R>=0`, `V_x=-sR<=0`, and `V<=1`. Consequently

\[
 0\le R(x,a)\le h/k_-,\qquad
 V(x,a)=1-\int_a^x s(t)R(t,a)\,dt
       \ge1-\frac{s_+h^2}{2k_-}
       \ge\frac{283}{608}>0.                                    \tag{5}
\]

This contradicts the first zero. Thus the positive-flux bootstrap is valid
throughout the interval. There was no prior assumption that the Green
function of a positive-frequency equation must stay positive arbitrarily
long; the explicit short-window margin is essential.

Dividing the now-positive lower flux by `k_+` yields

\[
 R_x(x,a)\ge\frac{1-s_+h^2/(2k_-)}{k_+}
      \ge\frac{1415}{12768}>\frac1{10}.                           \tag{6}
\]

The last strict margin is `691/63840`. Integrating the **pointwise**
positive lower flux, instead of replacing it by its endpoint minimum,
gives the sharper kernel bound

\[
 R(x,a)\ge\frac h{k_+}-\frac{s_+h^3}{6k_-k_+}
       \ge\frac{7495}{38304}h>\frac3{16}h\quad(h>0).              \tag{7}
\]

The final coefficient margin is `313/38304`; at `h=0` both kernel bounds
hold by equality and continuity. For later use, the same bootstrap also
supplies `R_x<=1/k_-` and `R<=h/k_-`.

## 3. Independent direct Volterra lower-bound route

For comparison, integrating (5) before estimating its two terms gives

\[
 R(x,a)=\int_a^x\frac{dt}{k(t)}
 -\int_a^x\frac{dt}{k(t)}\int_a^t s(z)R(z,a)\,dz.
\]

Applying `R<=h/k_-` directly to the second integral yields the weaker but
still sufficient inequality

\[
 R(x,a)\ge h\left(\frac1{k_+}-\frac{s_+h^2}{6k_-^2}\right)
      \ge\frac{46385}{242592}h>\frac3{16}h.                       \tag{8}
\]

Its last margin is `899/242592`. Both derivations are retained in the
exact rational replay. The stronger coefficient (7) uses the fact that
the lower flux was proved positive before dividing by `k_+`.

Increasing the allowed interval can invalidate (5). The executable
longer-window control makes that sufficient margin negative; it is not
relabelled as a proof that the actual longer-window Green function is
negative or dynamically unstable.

## 4. Smooth pulse class and retarded endpoint response

Use zero initial hidden data at the fixed left endpoint,

\[
 Q(-1/4)=0,\qquad Q_x(-1/4)=0.                                    \tag{9}
\]

Let `eta>0` and let `q` be any smooth nonnegative pulse with compact
support in `(-1/4,-1/8)`, `q<=eta`, and

\[
 \int q(a)\,da\ge\eta/16.                                        \tag{10}
\]

Variation of constants for (2), with the normalization (4), gives

\[
 Q(x)=\int_{-1/4}^x R(x,a)s(a)q(a)\,da,
 \qquad Q_x(x)=\int_{-1/4}^x R_x(x,a)s(a)q(a)\,da.                 \tag{11}
\]

There is no boundary term in the second formula because `R(x,x)=0`.
Every source point at final time has `-a>=1/8`. Using the conservative
bounds `R>=3h/16`, `R_x>=1/10`, `s>=42` and (10) proves

\[
 Q(0)\ge\frac{63}{1024}\eta,\qquad
 Q_x(0)\ge\frac{21}{80}\eta.                                     \tag{12}
\]

The sharper Green constants give still larger lower bounds; (12) is kept
as a simple common calibration. These inequalities hold for every finite
delta in (1), and for every pulse in the specified class.

The class is nonempty without a transcendental integration estimate.
Define `rho(t)=exp(-1/t)` for `t>0` and zero otherwise, and the smooth step
`S(t)=rho(t)/(rho(t)+rho(1-t))`. It is zero through `t=0` and one from
`t=1`; its denominator is everywhere positive. Choose

\[
 L=-31/128,\quad R_0=-17/128,\quad w=1/64,\qquad
 q(x)=\eta S((x-L)/w)S((R_0-x)/w).                                \tag{13}
\]

Its support closure is strictly inside `(-1/4,-1/8)`. Its unit-height
plateau is `[-29/128,-19/128]`, of length `5/64`, so its area is at least
`5eta/64>eta/16`. To verify smoothness, for positive `t` write
`D_t^n rho=exp(-1/t) P_n(1/t)`. The exact recurrence
`P_(n+1)(z)=z^2[P_n(z)-P_n'(z)]` gives a polynomial of degree `2n`.
Exponential decay dominates that polynomial at zero, proving that every
endpoint derivative vanishes. The quotient and product in (13) are
therefore globally smooth, including the plateau joins.

## 5. Small amplitudes, identical germs, and unselected inverses

The upper bounds following (7), together with `q<=eta`, give throughout
the whole window

\[
 0\le Q\le\frac{975}{2432}\eta,\qquad
 0\le Q_x\le\frac{325}{152}\eta.                                 \tag{14}
\]

For the first constant, use
`integral_(-1/4)^(-1/8) (-a) da=3/128`; monotonicity of the response gives
its maximum at the final time. For the second, the support interval has
length `1/8`. Thus an independent small amplitude makes both metric
anisotropies small in the linear problem uniformly in the admitted delta.
This is not a bound on a nonlinear remainder or on physical time derivatives.

The zero input and (13) have the same entire final `q` germ: each is
identically zero on an open neighborhood of `x=0`. They share the same
background coefficients and the same zero initial hidden data. The former
has zero retarded hidden response, whereas the latter obeys (12). No map
of that final germ alone can exactly reproduce both hidden phase outputs.
A map calibrated to the zero input has at least the stated pulse error;
an arbitrary common final output for the two histories incurs at least
half their separation on one of them. A stateful/nonlocal map is not
excluded by this elementary comparison.

The zero-data prescription is indispensable. With unrestricted homogeneous
hidden data, one can add another solution of the homogeneous equation and
change or cancel the endpoint output. This theorem concerns a declared
retarded problem, not all inverses or all histories allowed by another
matching contract.

## 6. Physical and matching scope

The full physical window lasts `tau*sqrt(delta)/4`. The derivative in (12)
means

\[
 Q_T(0)\ge\frac{21\eta}{80\tau\sqrt\delta}.                      \tag{15}
\]

For fixed eta the input and response vary on a shrinking physical time
scale. They are not a fixed low-frequency test band below a heavy gap.
Allowing eta to shrink preserves the normalized lower bounds while the
absolute metric amplitudes shrink; it does not supply a cutoff or an
independent nonlinear error estimate.

The input here is a **prescribed physical-metric history**, not an
externally conserved matter stress and not the internal free-chi source.
Maintaining `q=0` after the pulse while `Q` remains nonzero would generally
require driving in the physical-g equation. No on-shell g solution or
physical source-cost statement is inferred. The root-owned physical-g
Euler readout is a separately normalized diagnostic, not a replacement
for this distinction.

All bounds solve the exact finite-delta linear equation in divergence
form. They neither use the limiting own-f coefficient 16 as a finite-delta
substitute nor transplant the coupled-relative coefficient 80 into this
prescribed-q problem. No conclusion about every second-order EFT,
original CD/M1 matching, scalar/vector health, quantum effects, ultraviolet
completion or closure of P8 follows.
