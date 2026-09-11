# Actual weighted Hamiltonian and original forces

Write vhat for the regular clock-chart log-scale perturbation,
n for the physical lapse, and p=Pi/a0^3. The reference scale
is a0=(1+t^2)^2, H=4t/(1+t^2). S6.180's full-action derivation
retains its lapse-dependent canonical boundary before
differentiating. Its correct fixed-charge Routh density is

    L0=-3vhat'^2+6theta n vhat'+(J-3theta^2)n^2
       -L n vhat-(9/2)ell^2 vhat^2,
    delta=1/[2(1+t^2)^3], ell=1/[10(1+t^2)^6],
    w0=ell(3delta-1), L=3ell w0,
    theta=H-t/(1+t^2)^4.

The polynomial P=800(1+t^2)^18 J is stored explicitly in
hamiltonian.py and independently checked against the
literal full-target J. All its nonzero coefficients are
positive and even-degree, with constant1215.

The physical scale is v=vhat+delta n. Adding external
force density -gN n-gV v therefore gives transformed lapse
force gN+delta gV, not just gN. These are the original
gauge-fixed metric Euler forces; the clock equation is
the corresponding Ward equation, not an independently
chosen additional force. No M1 force is added.

The momentum is p=-6vhat'+6theta n. Its Hamiltonian,
including the external force, is

    Hred_before=-p^2/12+n(theta p+L vhat)-J n^2
               +(9/2)ell^2 vhat^2+(gN+delta gV)n+gV vhat.

Its algebraic constraint is

    n=(theta p+L vhat+gN+delta gV)/(2J).

Only J is inverted. The velocity-chart coefficient
J-3theta^2 need not stay away from zero and is never a
denominator. Eliminating n gives

    Hred=-p^2/12+(theta p+L vhat+gN+delta gV)^2/(4J)
         +(9/2)ell^2 vhat^2+gV vhat.

Since p is normalized by the time-dependent a0^3, Hamilton's
equations are vhat'=partial_p Hred and
p'=-partial_vhat Hred-3Hp. Their explicit matrices A,B
are stored and checked both from this Legendre transform
and the original two Euler equations with sources.

The original M1 equation reconstructs
sigma'=-w0 n-3ell vhat. Its independent charge variation is
zero; the background charge is not removed. The normalized
background obeys ell'+3H ell=0. Prepared zero data for the
metric and sigma give the unique smooth reference response.
This first-order ODE and algebraic readout are regular
through the bounce. No particle-ghost or cone conclusion
is inferred from a sign in this reduced homogeneous chart.
