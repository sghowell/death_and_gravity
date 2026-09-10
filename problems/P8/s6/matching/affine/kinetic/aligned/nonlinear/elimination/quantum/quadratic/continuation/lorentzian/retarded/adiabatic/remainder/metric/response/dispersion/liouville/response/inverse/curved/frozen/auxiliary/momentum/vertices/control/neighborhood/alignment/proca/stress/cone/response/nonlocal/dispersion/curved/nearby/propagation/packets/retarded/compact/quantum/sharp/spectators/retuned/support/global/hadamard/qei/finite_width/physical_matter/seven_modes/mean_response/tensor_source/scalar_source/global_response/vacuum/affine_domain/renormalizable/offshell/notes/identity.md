# Full four-dimensional pointwise identity

All derivatives use the flat +--- metric. Set

    E=Box(phi)+phi, X=(dphi)^2,
    Z=dphi.Hess(phi).dphi,
    P=phi E, J=phi^2, K=Box+2,
    c=G^2/D^2=4lambda^2/gamma.

The exact constant and linear centered-resolvent terms are

    L01=c J^2/12-c J KJ/8.

The next terms have coefficients lambda/4 and -gamma/8.
Thus L4_trunc=L01+(lambda/4) J K^2J-(gamma/8) J K^3J.
These coefficients are checked from the literal polynomial
potential and full heavy inverse, including lambda4.

An explicit cubic field polynomial is

    R=-c phi^3/6+lambda(2phi X+phi^2 E)
      +gamma(2Z+phi X+phi^3/3)
      -gamma phi KX-(gamma/2)phi K(phi E).

Using contravariant gradients, take the current

    j=-c phi^3 grad(phi)/12
      +(lambda/4)[J grad(KJ)-KJ grad(J)]
      -(gamma/8)[J grad(K^2J)-K^2J grad(J)
                 +4(X grad(P)-P grad(X))]
      -(gamma/2)X grad(X)
      -gamma(phi X+phi^3/3)grad(phi).

Then the literal identity is

    L4_trunc-L4_target-E R-div(j)=0.

The written derivation uses KJ=2(X+P), the two Green identities
J K^2J-(KJ)^2=div[J grad(KJ)-KJ grad(J)] and
J K^3J-(KJ)K^2J=div[J grad(K^2J)-K^2J grad(J)],
and the cross identity X KP-P KX=div[X grad(P)-P grad(X)].
Also (dX)^2=4L4. The remaining massive gamma terms satisfy

    2gamma L4-gamma X^2
      -[-2gamma(L3-L4)-gamma phi^4/3]
      =gamma E(2Z+phi X+phi^3/3)
       -gamma div[(phi X+phi^3/3)grad(phi)].

The native calculation reconstructs every derivative using
independent multi-index symbols in four dimensions, not a
homogeneous ansatz or random numerical sample. There are
210 symmetric jets through order six. It checks the full
pointwise polynomial and each displayed intermediate identity.

The target scalar quartic is independently read from S6.109's
actual local expansion at kappa=n/gamma. Its removable
A3=-2n, A4=2n and A5=0 at the vacuum give -2gamma(L3-L4)
after canonical normalization. The same mass-one quadratic
is checked. No new finite-M graviton identity is asserted.
