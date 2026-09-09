# Finite reference-state difference inequalities

Let W_ref be an S6.99 state, W_omega any state in the same
generalized Hadamard class, and D_a finitely many real smooth
differential operators on the ACTUAL observable O. Along a subclock
curve define K_omega=sum_a pullback[(D_a tensor D_a)W_omega].

The contraction test in notes/observers.md gives the pullback and
its separated time-frequency signs. It is positive type.
K_omega-K_ref is smooth and symmetric: it has the same commutator,
and its two possible wavefront signs are disjoint.

For real compactly supported g let
F_omega(alpha)=K_omega(g*exp(-i*alpha*s),g*exp(i*alpha*t)).
Positivity gives F_omega(alpha)>=0. Symmetry and Fourier inversion
give

    integral g^2 (K_omega-K_ref)(s,s) ds
       = (1/pi) integral_0^infinity [F_omega(alpha)-F_ref(alpha)] d_alpha
       >= -(1/pi) integral_0^infinity F_ref(alpha) d_alpha.

The reference integrand decreases faster than any power in this
opposite-sign Fourier direction, so the bound is finite.
This is the positive-type argument of
[Fewster, Section 4](https://arxiv.org/pdf/gr-qc/9910060), with
the actual subclock and generalized-state hypotheses supplied here;
no Klein--Gordon equation is substituted.

Weights such as one half are absorbed into the D_a. In particular
this applies to O itself, real directional derivatives, and the
four derivative squares defining the fixed-background test energy.
The exact finite-width bound is this reference functional. The
small-sampling coefficient computed separately does not replace it.
