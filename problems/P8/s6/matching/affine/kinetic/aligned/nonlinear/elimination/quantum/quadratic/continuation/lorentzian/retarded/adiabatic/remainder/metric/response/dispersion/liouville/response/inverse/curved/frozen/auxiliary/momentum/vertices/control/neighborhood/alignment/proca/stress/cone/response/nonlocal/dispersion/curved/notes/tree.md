# Literal tree and the missing classical direction

The S6.80 source removal and S6.81 ordinary-Proca replacement preserve
the original background and quadratic action; their independent full
quadratic bridge checks are pinned upstream. S6.82's fixed profile is
counted with the quantum response, not inserted into this classical
tree and then counted again.

Start with the actual normalized homogeneous tree density

    -3 vhat'^2+(J_e+w0^2/2-3Theta^2)n^2+6Theta n vhat'
    +s'^2/2+w0 n s'-3ell vhat' s.

Here w0 is the fixed background matter coefficient; it is not the new
metric source w. Since ell'=-3H ell, the matter equation is

    [a^3(s'+w0 n+3ell vhat)]'=0.

Prepared data select zero charge perturbation, hence
s'=-3ell vhat-w0 n. Substitute this in the actual metric equations,
or equivalently vary the fixed-charge effective density

    L_eff=-3 vhat'^2+6Theta n vhat'+(J_e-3Theta^2)n^2
          -3ell w0 n vhat-(9/2)ell^2 vhat^2.

The upstream literal action and both reduction routes are replayed
independently. Merely substituting an on-shell matter velocity into the
unreduced Lagrangian is not used as a variational rule.

Now n=eta' and vhat=w+H eta-delta eta', where delta=1/(2h).
The complete density is substituted before taking any Euler derivative.
For a density depending on eta'' use the weighted Euler operation

    dL/deta -(partial_u+3H)dL/deta'
      +(partial_u+3H)^2 dL/deta''.

Keep all derivatives of H,Theta,delta,ell,w0,J_e. A second route first
takes the actual physical two-current operator and then applies
E_eta=-(partial_u+3H)E_N+H E_Z, E_w=E_Z. The two routes agree exactly.

The complete local derivative inventory has orders [[4,3],[3,2]].
Its top entries are

    T_eta,eta: -6 delta^2 partial_u^4,
    T_eta,w:   +6 delta partial_u^3,
    T_w,eta:   -6 delta partial_u^3,
    T_w,w:     +6 partial_u^2.

The negative fourth coefficient is not interpreted as a healthy
particle kinetic term. It is the nonzero coefficient needed for this
prepared integral-equation reduction. The scalar degree-of-freedom
or stability interpretation is outside this inversion theorem.

On |u|<=1/2, delta=1/[2(1+u^2)^3] lies in [32/125,1/2].
Its derivative -3u/(1+u^2)^4 proves the continuous extrema, with the
endpoint and center values checked exactly. Thus for A=-6delta^2,

    |A|>=6144/15625, |A^-1|<=15625/6144.

No division by H,Theta,frequency,momentum or a vanishing bounce rate
appears. The overall classical L^2 has already been factored out of
T; gamma carries the relative quantum normalization.

The physical sources and matter are recovered by

    n=eta', v=w+H eta,
    s'=-3ell(w+H eta-delta eta')-w0 eta',

with zero past displacement. These formulas retain the original
negative reconstruction sign and determine the original matter
equation and charge condition, not an independently chosen solution.
