# Both principal charts and continuous margins

Keep q=k_com²/a² and physical cosmic u. All matrices below
are divided by the common positive a³. The units do not
change the physical matter cone.

## Unitary chart away from Theta=0

The lower-scalar addition contributes no q-weighted term after
the unchanged scalar-shift constraint. Thus its gradient matrix
is the old one, G=K_old, while

    K=K_old+diag(delta_J/Theta²,0),
    K_old=[ (J+w²/2)/Theta²   w/(2Theta) ]
          [ w/(2Theta)               1/2 ].

The positive-square form gives K>0 and G>0 for J>0,delta_J>0.
The difference K-G is positive semidefinite of rank one.
An exact determinant calculation gives

    det(G-c²K)/det(K)=(1-c²)[J/(J+delta_J)-c²].

The longitudinal, transverse and tensor principal blocks remain
those of the aligned model. Their physical front speeds are one;
a massive finite-q phase velocity is not used as a cone test.

## Regular gamma chart at the crossing

Start with the same regular Hamiltonian and canonical momenta
P_b=2a³q*v,P_s=a³*p_m. The canonical Hamiltonian includes
the moving-generator term -H*b*P_b. Write it as

    H=P^T A P/2+P^T B Q+Q^T C Q/2.

The velocity action has alpha=A^-1,beta=-A^-1 B,
gamma=B^T A^-1 B-C. Its Euler gradient is beta'-gamma, not
-gamma alone. At high physical q the exact rational calculation
gives

    K_gamma=[ (J+delta_J+w²/2)/Lambda²  -w/(2Lambda) ]
            [ -w/(2Lambda)                        1/2 ],
    beta/(a³q) -> diag(-2Theta/Lambda,0),
    gamma/(a³q) ->
        [ -2Theta(2H Lambda-Theta)/Lambda²   -ell ]
        [ -ell                               -1 ].

The leading beta is symmetric and gamma has no q² term.
Since (a³q)'=H*a³q, the principal gradient is

    G_gamma,11=(-Lambda*Theta'+H*Lambda*Theta
                 +Theta*Lambda'-Theta²)/Lambda²,
    G_gamma,12=ell/2, G_gamma,22=1/2.

The actual background identities, independently replayed, are

    w=-ell*Lambda,
    J+w²/2=Theta(H Lambda+Lambda')-Lambda*Theta'-Theta².

Thus G_gamma is the old positive K_gamma and the added kinetic
piece is diag(delta_J/Lambda²,0). Its characteristic polynomial
is the same as in the unitary chart.

These are genuine local principal limits on Lambda!=0:
the finite-q denominator is
D=q*Lambda²-J-delta_J-w²/2. On a compact neighborhood
where Lambda stays nonzero, its rational large-q expansion and
each needed time derivative have uniform remainder estimates.
The beta antisymmetric part has order one, so it does not add
an unaccounted order-q or order-sqrt(q) principal mixing.
This statement is about the classical characteristic symbol,
not an EFT validity estimate at arbitrarily high frequency.

At u=0, Lambda=-1/2,Theta=0,Theta'=3,H'=4. Differentiating
the full beta before taking its high-q limit independently gives

    G_gamma(0)=[6,1/20;1/20,1/2].

The finite-q positive velocity chart has q>6+16epsilon.
The regular phase Hamiltonian itself has no such pole. Omitting
the time derivative would falsely give G_11(0)=0.

## Chart coverage and all-time positivity

On |u|<=1/4, h<5/4 and Lambda=1-3/(2h)<=-1/5.
On 1/4<=|u|<=1/2,

    |Theta|=|u|(4-1/h)/(1+u²)>=3/5.

Hence these two regular principal charts cover the full compact
interval including the bounce. For any other finite time,
Theta is nonzero and the unitary chart applies. There is no
uncertified Theta=0 exception.

Let J*h²=P34(u)/[800(1+u²)^12]. Subtracting 1199/800
leaves a numerator with only positive even-power coefficients
and zero constant. Therefore J*h²>=1199/800 for all real u.
On |u|<=1/2, positivity of those coefficients and denominator
at least 800 give the explicit upper bound

    J*h² <=3817516721573/107374182400 <36.

Consequently, for epsilon>0,

    c_clock² >=1199/(1199+3200epsilon)>0,
    1-c_clock² >=4epsilon/(3817516721573/107374182400+4epsilon)
        on |u|<=1/2,
    fractional added light kinetic form <=3200epsilon/1199.

The last inequality follows directly from the old positive
square form and is independent of principal chart. For
epsilon=10^-6 it is below three parts per million, while
the compact speed margin exceeds 10^-7.

The free-matter direction, and the tensor/vector blocks, remain
luminal. This is not a strictly positive margin in every
direction and cannot absorb arbitrary unsigned quantum
principal perturbations. A negative epsilon, while small
enough to keep K positive, makes the clock superluminal;
that sign is an explicit negative control, not an allowed model.
