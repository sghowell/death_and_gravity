# Physical lapse recovery and the explicit local constant

A C0 bound for eta,w alone would not control n=eta'.
The actual first integrated row has the stronger local structure

    A(t)eta(t)+integral V_eta(t,s)z(s)ds=I4 g_eta(t).

Differentiate ONCE in output time without differentiating z(s).
With a0=sup|A^-1|,

    ||eta'|| <= a0 [||I3 g_eta||+C1 ||z||],

where C1 bounds |A'|, the row norm of V_eta(t,t), and the
integral of its first output-derivative norm.

The full current tree's derivative coefficients are exact
rational functions. Their derivatives through the orders needed
by the local primitive kernels, and A' through A'''', have
continuous rational coefficient enclosures on I.
All24 reconstructions are checked. The variable-fourth term
contributes -4A' to the kernel diagonal; no such contact is
dropped.

The two normalized quantum time-row coefficients are
-3H' P_physical/kappa (order zero in eta) and
-3P_physical/kappa (order one in w).
S6.176 gives |P|/kappa,|P'|/kappa<1e-770 and |H'|<=4.
The primitive kernels' output derivatives are bounded by an
additional (21/2)*1e-770; these terms have no diagonal contact.
The exact assembled rational bound is

    C1 <15000

(approximately14976.3097083, only an illustrative display).
This is a quantitative LOCAL reconstruction constant, not the
full inverse's C,K, exponential weight or stability bound.

For original physical force g=(g_N,g_v), the adapted force is

    g_eta=-(partial_t+3H)g_N+H g_v, g_w=g_v.

Prepared integration by parts removes the derivative of g_N.
With |H|<=8/5 and a unit-length slab,

    ||I4 g_adapted|| <=13/30 ||g||,
    ||I3 g_eta|| <=47/30 ||g||,
    ||v||=||w+H eta|| <=13/5 ||z||.

Thus continuous original forcing has a finite C0 inverse bound
for BOTH original metric sources. Smooth prepared forcing gives
smooth prepared solutions: shift derivatives to the prepared
input and to diagonal kernel derivatives, which retain integrable
majorants, and induct using the same bounded Volterra inverse.

Differentiating the integrated equations four times recovers the
retained response equations in the zero-past distribution algebra.
With the scale residual zero, the missing physical lapse residual
obeys (partial_t+3H)E_N=0; its prepared data force it to zero.
The source transformation loses no original equation.

Finally reconstruct the original matter displacement by

    sigma'=-3ell(w+H eta-delta eta')-w0 eta'.

Its zero initial displacement and zero charge perturbation are
retained; no matter datum is reset. Causality and uniqueness hold
on every initial subinterval. They do not exclude a large transient,
growing modes or frequencies outside an unknown EFT cutoff.
