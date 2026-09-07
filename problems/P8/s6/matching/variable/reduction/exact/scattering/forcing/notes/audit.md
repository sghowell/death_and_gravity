# Independent S6.36 source, phase and original-field audit

This audit derives the loading and readout independently and compares the
candidate to the immutable S6.35 `time_transfer`, `plane_wave_frame` and
`endpoint_map` APIs. The pinned report is
`c77d6d3e41d81c5d2bda33aeba201fe438620d872f3e9f2d49fde5469194641e`.
The tests include the distinct Fraction reconstruction in `audit.py` but
do not use its returned constants as their only evidence.

## 1. Actual source and a strictly positive limiting load

Let `L=1/100`, `a=L/2` and `0<delta<=10^-6`. On the fixed punctured
interval `[-L,-a]`, the parent branch gives `19/10<b,N<21/10`. Therefore

    1 < (19/10)^3/(21/10) < k=b^3/N
      < (21/10)^3/(19/10) < 5 < 25/4.

Put `m=A_profile/D=2 beta1_bar N/b^2` and `Y=sqrt(k) Q`.
Direct expansion of the own-metric equation gives

    (k Q_u)_u+k m(Q-q)=0
    <=> Y_uu+[m-(sqrt(k))_uu/sqrt(k)]Y=sqrt(k) m q.

Thus the canonical source is not just `q`. The parent remainder and pump
bounds give, throughout this interval,

    0 < 16/(delta_max+8L^2)-44 <= V <= 8/L^2+44 =80044,
    m >=16/(delta_max+8L^2)-44-22 >1/L^2.

For a nonnegative prescribed pulse, `sqrt(k)m q>=q/L^2`; it is harmless
that equality holds where the pulse is zero. Zero physical data at `-L`
are also zero canonical data, including the time-dependent normalization
jet. There is no incoming homogeneous field in this construction.

For the unit retarded Green function starting at a source time `s`, set
`h=u-s`. If `0<=V<=M`, the first-zero argument gives `0<G<=h` until a
possible first zero of `G_u`. On that interval its Volterra identity yields

    G_u>=1-M h^2/2,     G>=h-M h^3/6.

A strictly positive first bound excludes that first zero, closing the
bootstrap. At `M=80044`, `h<=L/4`, these bounds are

    G_u>=59989/80000>2/3,
    G/h>=219989/240000>9/10.

This `L/4` lag comes from the actual support, not the full initial-to-load
length `L/2`. The latter would fail this particular bootstrap. The explicit
smooth pulse has support `[-95L/128,-81L/128]`, strictly inside
`(-3L/4,-5L/8)`, and plateau width `5L/64>L/16`. Consequently its area is
greater than `eta L/16`, and every contribution at `u=-a` has lag at least
`L/8`. In particular the following *fixed stronger* lower bounds hold:

    Y_delta(-a)/eta >= (219989/240000)/128
                       =9/1280+3989/30720000,
    Y_delta,u(-a)/eta >= (59989/80000)/(16L)>1/(24L).

The proof of a strictly positive limiting load uses these uniform gaps,
not the invalid rule that a strict inequality remains strict under limits.
The denominator on this punctured interval obeys
`D>=2-2/(1+a^2)^4>0`. The parent's analytic implicit branch hence gives
ordinary smooth coefficient convergence as `delta` decreases to zero
here. Continuous dependence for the linear forced ODE supplies a fixed
limiting load. No coefficient convergence at the singular center is used.

## 2. The correct real phase and exact parameter sequence

Write the parent Jost coefficient as `A_J=ar+i ai` and `B_J=i beta`.
For `theta=rho asinh(sqrt(8)a/sqrt(delta))`, both endpoint wave frames
must be retained. With `V0=2^-1/2[[1,1],[i,-i]]`, the reference Cauchy map is

    V0 diag(e^(i theta),e^(-i theta)) T_t
       diag(e^(i theta),e^(-i theta)) V0^*,

not `T_t` alone and not the oppositely oriented radial-power transfer.
For `c=cos(2theta)`, `s=sin(2theta)`, its real matrix is

    [[ar c+ai s, ar s-ai c-beta],
     [-ar s+ai c-beta, ar c+ai s]].

The determinant is `|A_J|^2-beta^2=1`. For the nonzero real loaded limit
`w=(w1,w2)`, the first component is exactly

    c C+s S-beta w2,
    C=ar w1-ai w2,  S=ai w1+ar w2,
    C+i S=A_J(w1+i w2), C^2+S^2=|A_J|^2 ||w||^2.

Choose `theta0=arg[A_J(w1+i w2)]/2` modulo `pi` in `[0,pi)`.
Then `theta0+pi n` and `theta0+pi n+pi/2` give opposite scalar extrema.
The factor `1/2` in this offset is essential. The constant `beta` term
cancels between the pair. It is the oscillatory `A_J` part, not merely
nonzero `B_J`, that carries this phase dependence.

The exact positive parameter sequences are

    delta_n^+ =8a^2/sinh^2[(theta0+pi n)/rho],
    delta_n^- =8a^2/sinh^2[(theta0+pi n+pi/2)/rho].

They invert the actual asinh clock, rather than an asymptotic logarithm.
For `n>=2`, `pi>3`, `rho<4/3` and nonnegative offset give the argument
greater than `9/2`. Since
`sinh(9/2)>9/2+(9/2)^3/6=315/16`, both deltas lie below `10^-6`.
Their monotonic decrease to zero follows directly from the derivative of
`8a^2 csch^2(theta/rho)`. The offset depends on this one fixed source's
limiting load; it does not change the source or incoming data with delta.

## 3. Two finite errors, no unproved vanishing remainder

The smaller central half-width `a=1/200` requires a fresh budget. The
canonical perturbation integral is at most `33/25000`, the two removed
reference tails at most `3/800`, and division by `rho>5/4` gives

    J <507/125000<1/200.

The parent's source-free comparison proof then bounds the coefficient-map
error by `epsilon=4/399`. Unitary endpoint rotations preserve this norm.
This estimate is fixed, not `o(1)` at fixed physical width. The two phase
sequences can therefore contribute two different errors, each of size
`epsilon ||w||` in the limit. Their total loss is `2epsilon ||w||`, not
one error and not zero. This leaves

    2(|A_J|-epsilon)||w|| > (790/399)||w||.

For the aligned first component, one sequence is bounded below by
`(|A_J|-epsilon)||w||-beta w2+o(1)` and the other above by
`-(|A_J|-epsilon)||w||-beta w2+o(1)`. The `o(1)` here comes only from the
proved convergence of the punctured loading vector and endpoint factors;
it is **not** an assumed limit of the central coefficient-map error.
These one-sided bounds imply the required global limsup/liminf separation
even if neither actual output sequence converges. A bounded alternating
error is an explicit norm-admissible countercontrol to replacing the
comparison estimate by a vanishing error or claiming actual sequence limits.

## 4. Read out the original own-metric amplitude

The frozen exact endpoint map has first row
`Z1=sqrt(k/r_delta) Q`, where `r_delta=sqrt(a^2+delta/8)`. Thus

    Q(+a)=sqrt(r_delta/k_delta(+a)) Z_out,1.

It is the unchanged original `Q`, not an adjustable projection of a
two-component field. The full endpoint map retains `k_u`; its input
derivative is `Q_u=tau Q_T`. Its inverse first row nevertheless contains
no velocity coefficient. At the punctured limit, `w1=Y_*(-a)/sqrt(a)>0`.
The physical limiting readout factor therefore satisfies

    sqrt(a/k_*(+a)) ||w|| >=Y_*(-a)/sqrt(k_*(+a))
                            >(9eta/1280)(2/5).

Combining this with both transfer errors gives

    limsup Q_delta(+a)-liminf Q_delta(+a)
      > eta*(790/399)*(9/1280)*(2/5)
      = eta*237/42560 >eta/200.

The pulse and observation times are fixed in `u=T/tau`, with fixed tau.
This is a zero-data response to a prescribed physical-metric tensor
history on the off-shell own-f branch. It is not automatically a conserved
matter source or a solution of the physical-g field equation. The pulse is
smooth and compact, not an exact temporal bandlimit, and no uniform
low-frequency/instantaneous-stiffness hierarchy, EFT cutoff, vacuum
production, or UV conclusion follows from this response theorem.
