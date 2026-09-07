# Independent physical-tensor obstruction

This is a test of the **formal retained action**, not a high-frequency
test of the complete two-metric theory and not a cutoff calculation.
The metric remains the physical g to which the free scalar chi couples.
No equations of motion, metric redefinitions or homogeneous-mode
restrictions are used to identify different actions.

## 1. A literal off-shell metric calculation

Take an arbitrary homogeneous clock phi(t) and

    g = diag(1, -exp(q(t)), -exp(-q(t)), -1).

This is a trace-free tensor variation about q=0 with exactly unit volume.
It need not solve any background equation. The coordinate definition of
the connection and curvature in the frozen P8 formulation gives

    R_00 = -q_t^2/2,
    R_11 = exp(q) q_tt/2,
    R_22 = -exp(-q) q_tt/2,
    R_33 = 0,
    R_B = -q_t^2/2.

Thus, without expanding or integrating by parts,

    R_mn R^mn - R_B^2/3 = q_tt^2/2 + q_t^4/6.

The root-authored `tensor.py` computes the full connection and Ricci
tensor from the metric. It does not import either stationary-metric
reduction or scalar coefficient dictionary code.

For p=phi_t and z=phi_tt the five literal clock contractions are

    L1 = z^2 + p^2 q_t^2/2,
    L2 = z^2,
    L3 = L4 = p^2 z^2,
    L5 = p^4 z^2,
    X = p^2,  box(phi) = z.

Consequently **every** action in the frozen quadratic scalar-tensor
basis, with arbitrary coefficient functions, has tensor density

    L_T = (-F2+X A1) q_t^2/2 = G_T q_t^2/4.

Neither the Ia relations nor the clock or metric field equations were
needed. The coefficients are functions of t on this off-shell test
because the clock and X do not depend on q. Its tensor Euler equation
is second order even when those coefficients vary.

## 2. The retained higher operator cannot be a boundary

The stationary-f calculation retains

    kappa(phi) (R_mn R^mn - R_B^2/3),
    kappa = M^4 r^3/(4 beta1).

Its tensor second variation contains kappa q_tt^2/2. The Euler variation
is

    kappa q_tttt + 2 kappa_t q_ttt + kappa_tt q_tt.

All time derivatives of kappa are retained. The leading fourth-order
coefficient is kappa, nonzero on the admitted positive branch. At the
center it is M^2 tau^2 c(c-2)/16 for 2<c<=4. The endpoint c=2 is not an
action in this result.

The clock-curvature cross terms in the full four-derivative reduction
have at most a second derivative of q multiplied by first derivatives;
the explicit scalar dictionary reduces them to the above second-order
tensor form. They cannot cancel the pure fourth-order coefficient.
Any total divergence has identically zero Euler variation, so equality
to a quadratic-DHOST action modulo boundaries is impossible in this
unchanged physical metric. A scalar-only invertible clock redefinition
does not introduce metric second derivatives into its Hessian and cannot
change this conclusion.

As an additional boundary control, the compact H0^2 variation
q=(1-t^2)^2 on [-1,1], zero outside, has q=q_t=0 at both endpoints and

    integral q_tt^2/2 dt = 64/5 > 0.

Smooth compact approximations preserve a positive value. This control
is not substituted for the pointwise Euler-symbol argument.

## 3. Independent source and matter checks

For constant proportional ratio, let the physical and hidden Einstein
weights be G and F, with positive relative spring nu. Keeping the g
source attached to g, the hidden field's own equation gives the exact
Schur kernel

    (G+F)D - F^2 D^2/nu + F^3 D^3/[nu(nu+F D)].

Here F=M^2 r^2 and nu=2 beta1 r. Hence the retained curvature coefficient
is kappa=F^2/(2nu), and the Weyl-square coefficient is kappa/2. The
physical g response is

    1/[(G+F)D] + F/[G(G+F)(D+m_FP^2)],
    m_FP^2=nu(G+F)/(GF).

The test derives both expressions from the literal two-by-two kernel.
Solving the other metric equation or dropping the physical source
weight does not give this dictionary. This constant-coefficient control
does not establish a rolling expansion window.

For the matter action sqrt(g) Y/2, an inverse-metric variation Z gives

    delta S_chi = (1/2) integral sqrt(g) T_chi,mn Z^mn,
    T_chi,mn = chi_m chi_n - g_mn Y/2.

A separately implemented rank-one metric variation verifies the sign
from the exact determinant lemma. Rank-one symmetric matrices span the
metric variations, so this is also a polynomial test of the full first
variation. The explicit curvature-eliminating inverse-metric change in
the companion dictionary therefore generates

    (kappa/Q) [R_mn chi^m chi^n - R_B Y/3],  Q=M^2(1+r^2).

It is not the free M1 matter action. If leading Einstein equations are
used as another operator change, its pure-chi contact is
2 kappa Y^2/(3Q^2); mixed clock contacts must also remain. We do not use
this on-shell rewriting as permission to drop source contact terms.

## 4. Scope

The fourth-order equation belongs to a **finite formal truncation**.
It is not a diagnosis of an extra ghost in the untruncated parent, nor
proof that this operator is important at every admitted physical
frequency. The adopted S6 contract permits higher operators with a
controlled error budget. This result establishes exact action
inequivalence and a necessary matching obligation. The companion
clock-normalized coefficient defect states the explicit norm in which
a small omitted remainder cannot repair the target. Other norms,
resummations, state restrictions or altered matter frames require
their own maps and errors; none is silently excluded or accepted.
