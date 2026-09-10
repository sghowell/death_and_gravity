# Actual derivative-map composite mixing

Use +--- signature, mass one, E=(Box+1)phi, X=(d phi)^2,
Z=d phi Hessian(phi) d phi and K=Box+2. The literal S6.111 map is

    R = -c phi^3/6 + lambda(2 phi X + phi^2 E)
        + gamma(2Z + phi X + phi^3/3)
        - gamma phi K(X) - gamma phi K(phi E)/2.

Let T_d be the coincident scalar Gaussian tadpole. In dimensional
regularization, polynomial loop integrals without a denominator
vanish. Dividing p^2 and p^4 by p^2-1 therefore gives I2=I4=T_d.
Odd moments vanish. The rank-two and rank-four moments are fixed
by Lorentz covariance and their traces. In particular,

    <d_mu phi d_nu phi> = eta_mu_nu T_d/d,
    <phi d_mu d_nu phi> = -eta_mu_nu T_d/d.

The rank-four tensor is the sum of three metric pairings times
T_d/[d(d+2)], with the derivative signs retained. Contractions of E,
or derivatives of E, against a field derivative reduce to scaleless
contacts and vanish in this prescription.

Let W denote one Wick pair contraction, leaving one field. Product
rules and these tensors give

    W[phi^3]/T_d       = 3 phi,
    W[phi X]/T_d       = phi,
    W[phi^2 E]/T_d     = E,
    W[Z]/T_d           = Box phi/d,
    W[phi K(X)]/T_d    = 2 phi - 4 Box phi/d,
    W[phi K(phi E)]/T_d = (Box+1)^2 phi.

For the fifth line, the two internal contractions in Box(X)
cancel: 2< Hessian(phi)^2 > + 2<d phi d Box(phi)>=0.
The two cross contractions with the outer phi yield -4 Box(phi)/d.
For the last line, expanding K(phi E) gives Box(E)-E+2E
after contraction, hence (Box+1)E. These derivations hold in
continued d before any pole subtraction.

Independently, native code expands the actual four-dimensional
polynomial in S6.111's 210 symmetric jets. It contracts each pair
in every cubic monomial. For multiindices a,b of total rank 2r,
the covariance divided by T_4 is

    (-1)^(r+|b|) sum(metric pairings)/[4*6*...*(4+2r-2)].

With diagonal metric, only even multiplicities per coordinate
survive. The pairing count is the product of the odd double
factorials. The result agrees with all six expressions above and
with the complete literal R, not merely a field-degree surrogate.

Writing B for the free Box eigenvalue, the assembled mixing is

    W[R] = T_d a_d(B) phi,
    a_d(B) = -c/2 + lambda(B+3) + 6 gamma B/d
             - gamma(B+1)^2/2.

Set d=4-2epsilon and use MSbar reference scale one. The mass-one
tadpole, including its continued dimensional factor, is

    T_d = exp(EulerGamma epsilon) Gamma(epsilon-1)/(16 pi^2)
        = -(1/epsilon+1)/(16 pi^2) + O(epsilon).

Native Gamma expansion checks both terms. Since
a_d=a_4+(3 gamma B/4)epsilon+O(epsilon^2), subtracting the pole of
the product gives

    W_MS(B) = -[a_4(B)+3 gamma B/4]/(16 pi^2).

At B=-1,

    w = W_MS(-1) = (c/2-2 lambda+9 gamma/4)/(16 pi^2).

Setting d=4 before subtracting would incorrectly omit
3 gamma/(64 pi^2) from this on-shell value. A negative control
retains that nonzero difference. This defines an explicit linear
composite-operator subtraction, not a transfer of any rolling-state
or finite-gravity prescription.

The inverse field is Psi=Phi-R(Phi)+O(Phi^5). The omitted inverse
term needs at least two Wick pairs to leave one particle. Adding
a parent quartic vertex, or two parent cubic vertices, to the
R one-particle graph also needs at least two loops. Thus the free
contraction is the complete one-loop linear mixing in this scheme.
Odd light parity preserves the zero vacuum expectation value.
Additional off-shell composite counterterms may be needed for
other correlators; their existence is not replaced by this
one-particle calculation.
