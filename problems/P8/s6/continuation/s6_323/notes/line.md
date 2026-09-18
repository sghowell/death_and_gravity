## Definition and universal identity

Let q_i=w_i n_i be future-null real rays, with w_i>0 and W<=1/8.
Each h_i is complex physical spatial TT with Frobenius norm<=1.
Use all-outgoing signed p, p^2=1. Strip the common canonical factor
kappa^(-|S|/2), and set j_empty=1. The exact scalar-line recursion is

 j_S = sum_nonempty_B_subset_S [
       -V_B(p+Q_(S minus B),-p-Q_S)*j_(S minus B)] / D_S,
 D_S=(p+Q_S)^2-1.

It has 1,3,13 trees at multiplicities1,2,3. No pure-soft graviton
branch is included in this line current. Every seagull and Phi2-h^3
contact is included. Define its connected normalized third coefficient

 U=abc*(j_123-j_12*j_3-j_13*j_2-j_23*j_1+2*j_1*j_2*j_3).

For fixed p this also equals the anchored threefold rectangle of
abc*j_123: its i-face is (w_j*w_k*j_jk)*(w_i*j_i), and every
pairwise-face intersection is the same product of three singleton
factors. This face identification follows directly by isolating the
terms with a singleton-i scalar cut; all others have no 1/w_i pole.
The compatible simultaneous faces additionally follow from the bound
below. This identification is NOT automatically asserted after
allowing p to depend on the three radiation energies.

Write d_i=2p.n_i and z_ij=2n_i.n_j. The seven original denominators are

 D_i=w_i*d_i,
 D_ij=w_i*d_i+w_j*d_j+w_i*w_j*z_ij,
 D_123=sum(w_i*d_i)+sum_pairs(w_i*w_j*z_ij).

The scalar one-graviton vertex, after using only the physical TT
conditions of its own external graviton, has the form

 L_i + sum_j_in_previous(m_ij*w_j+n_ij*w_j^2)
     + t_i*w_j*w_k  [if two previous leaves].

The two-graviton vertex has a constant, terms linear in its two own
energies, and, when a previous third leaf k is present, terms
w_k, w_k*w_i, w_k*w_j, w_k^2. The three-graviton vertex at base p is
constant plus linear terms in the three energies.

Treat ALL49 coefficients in these vertex/cut forms as independent
polynomial variables. Build the13 line terms, subtract the9 pair-times-
single terms, and add twice the triple-single product. Cross-multiply
the product of the seven ORIGINAL scalar cuts. The numerator P has
62 energy monomials; each is divisible byabc; its minimum degree is6
and maximum degree10. Thus P=abc*Q with every Q monomial of degree>=3.

This is an exact identity over a polynomial coefficient ring, not an
angular fit, numerical sampling argument, or altered GCD computation.
After normalizing byabc the denominator is

 d_1*d_2*d_3*D_12*D_13*D_23*D_123.

## Literal original-action specialization

The frozen action gives
 V(p+k,-p-k-q;A)=-2(p+k).A.(p+k+q)
                  +tr(eta*A)*[(p+k)^2-1+(p+k).q].
For physical TT A, tr(eta*A)=0 and A*q=0, so -V=2(p+k).A.(p+k).
Consequently L_i=2p.A_i.p, m_ij=4p.A_i.n_j,
n_ij=2n_j.A_i.n_j and t_i=4n_j.A_i.n_k.

For two or three fields put
 M=eta*(literal density-inverse coefficient)*eta
 and d=(literal determinant-density coefficient).
Their stripped vertex is d-(p+k).M.(p+k+Q_B).
Expanding this bilinear form gives exactly the generic degree pattern.
No TT, transversality or Ward premise is made for this density matrix.

The generic one-field identity and arbitrary-symmetric-M two/three
field coefficient identities were checked independently. At16 exact
fixtures (eight real/complex physical polarization combinations for
each sign of p), the complete original13-tree cumulant times its
original seven-cut denominator equals the specialized generic P.
These fixtures include three noncollinear rational directions, not
only the earlier orthogonal example.

## Uniform real bound

For either sign of p, the Doppler factor obeys |d_i|>1/2 and |d_i|<8.
Indeed E<=2 gives E-|p_space|>=2-sqrt(3)>1/4.
For any nonempty radiation subset with energy v<=1/8,

 |D_S| >= 2*(E-|p_space|)*v-v^2 > 3v/8.

The minus v^2 is only needed for the past-directed incoming p; for
future-directed p it is an enlargement of the bound. There is no
internal massless angular denominator on this scalar line.
Also |z_ij|<=4. Energy-polynomial coefficient-l1 bounds for D_i,
D_ij,D_123 are therefore8,20,36, respectively.

Each scalar momentum component has energy-polynomial coefficient-l1
at most2+3=5. With unit fields, the frozen literal density bounds are
K_r=(6,96,2400) and d_r=(4,48,960). Entrywise convolution therefore
bounds each vertex coefficient-l1 by16*25*K_r+d_r, namely

 V1=2404, V2=38448, V3=960960.

This is a coefficient-l1 estimate, not an inference from a real-point
supremum. The two scalar momentum factors contribute16 contractions
of coefficient-l1 at most25. Metric coefficients are energy-independent.

There are23 signed summands with total absolute weight24. Each
vertex product is bounded by max(V1^3,V1*V2,V3). Each missing-cut
product is bounded conservatively by8^3*20^3*36. Thus

 ||P||_coef <= 24*8^3*20^3*36*max(V1^3,V1*V2,V3) = C_P.

Removing its commonabc does not change its coefficient-l1. Since
every Q monomial has degree>=3 and0<W<=1/8<1, |P|<=C_P*abc*W^3.
The normalized denominator has modulus at least

 (1/2)^3*(3/8)^4*(a+b)*(a+c)*(b+c)*W.

Using
 (a+b)*(a+c)+(a+b)*(b+c)+(a+c)*(b+c)=W^2+ab+ac+bc>=W^2,
one obtains, restoring the canonical factor,

 |U_physical| < 1e24/kappa^(3/2)
                  *abc*[1/(a+b)+1/(a+c)+1/(b+c)].

The exact conservative coefficient is59670991094513926144000/3<1e24.
The bound gives compatible zero limits on every soft face and all
their intersections, uniformly in directions and the allowed p.
It is degree two in simultaneous soft scaling, with hierarchical
behavior retained. It is not a statement that U is a quadratic
polynomial or that its derivatives have a global-W Cauchy tube.

## Generic scalar pair split and parity

For pair energies a,b let d_i=2p.n_i,z=2n_i.n_j,
l=a*d_i+b*d_j and D=l+ab*z. From the literal scalar vertices,
with their coefficients defined in the connected-line proof,

 t_ij=ab*(R0+a*Ra+b*Rb)/(d_i*d_j*D),
 R0=-L_i*L_j*z+L_i*m_ji*d_j+L_j*m_ij*d_i+B_ij*d_i*d_j,
 Ra=L_i*n_ji*d_j+B_ij,qi*d_i*d_j,
 Rb=L_j*n_ij*d_i+B_ij,qj*d_i*d_j.

Define L_ij=ab*R0/(d_i*d_j*l). The exact remainder is

 R_ij=ab*[(a*Ra+b*Rb)*l-R0*ab*z]/(d_i*d_j*l*D).

Under p->-p, d_i,m_ij and B_ij,qi are odd; L_i,n_ij,B_ij,z
are even. Hence R0 is even, Ra/Rb odd and L_ij is odd.
The generic coefficient-ring calculation checks the entire identity
and this parity, without selecting angular data.

On the compact complex massive-momentum tube, and the pair-energy
tube of radii eta*(a+b), use
|L_i|<8, |m_ij|<7, |n_ij|<=2, |d_i|<10, |d_i|>1/4,
|B_ij|<14000, |B_ij,qi|<5000, |z|<=4.
The two B bounds follow directly from the literal density-inverse
coefficient96 and determinant coefficient48 with momentum
components<3:16*9*96+48<14000 and16*3*96<5000.

Writing u=a+b for the real center, |l|>u/8 and |D|>u/4.
The exact coefficient arithmetic gives
|R0|<1500000, |Ra|,|Rb|<600000, so

 |L_ij| < 2e8*|a'b'|/u,
 |R_ij| < 1e11*|a'b'|.

For fixed positive real a,b and complex p in a fixed spatial
polydisc of radius1e-6, the same estimates hold with strict margins.
The on-shell p0 branch remains near its real value; all Doppler
gaps stay>1/8, all scalar cuts>u/4. Applying Cauchy in p_space,
and keeping a,b real, gives a conservative
 ||grad_p(L_ij/u)|| < 1e16*ab/u^2.
The norm bound is below2e9*ab/u^2.
For the slightly larger complex a,b tube one may use the coarser
uniform2e9 and1e16 envelopes. These are pointwise p-derivative
bounds of an analytic massive branch, not a real-only gradient claim.

Consequently the real complete pair obeys
 |t_ij| < 2e10*ab/(a+b).
All these currents vanish on either pair face, with compatible
zero joint origin. Canonical powers are restored only after combining
the three emitted leaves.

## A real source pole excludes the wrong complex domain

For p=(1,0,0,0), q1=a*(1,1,0,0), q2=b*(1,0,1,0),
A=diag(0,0,1/2,-1/2), B=diag(0,1/2,0,-1/2), both tensors are
unit-bound physical TT. The literal connected pair has a nonzero
numerator at the scalar-cut root a=-b/(1+b). Take b=1e-30 and a
positive center a=b, with spectator c=1/16. The root lies strictly
within eta*W of the center, eta=1e-13. The production negative
control checks the nonremovable pole, not merely a zero denominator.

Consequently connected pair/triple currents do not receive a global-W
tube. The proof below uses real pair energies with c-only continuation,
or the relative-pair tube with strict scalar-cut margins. The hard-only
grouped kernels in notes/grouping.md are a different analytic object.

