# Smooth reference preparation on the original global clock

Define f(s)=exp(-1/s) for s>0, extended by zero for s<=0, and

    chi(s)=f(s)/(f(s)+f(1-s)),
    b(u)=1+chi(u+2)[a(u)-1],  a(u)=(1+u^2)^2.

The denominator never vanishes. The standard flat exponential
extension gives chi smooth on R, zero for s<=0, one for s>=1,
and 0<=chi<=1. Thus b>=1, b=1 on u<=-2, and b=a on u>=-1.
Both reference and physical metrics are smooth for every finite u.

Here is an explicit derivative bound, including the endpoint
regions rather than a numerical sample. On 0<s<1,

    chi'=chi(1-chi)[s^-2+(1-s)^-2].

On [1/4,3/4], each inverse square is at most 16, giving chi'<=8.
For 0<s<=1/4 put y=1/s>=4. Then
chi(1-chi)<=exp(4/3-y), (1-s)^-2<=16/9, and
y^2 exp(-y)<=16 exp(-4). Hence

    chi' <= (160/9) exp(-8/3)
          <= (160/9)/(1+8/3)=160/33<5.

Reflection covers the other endpoint. Therefore |chi'|<=8 globally.

On the transition strip [-2,-1], a<=25, |a'|<=40 and a-1<=24.
Consequently 1<=b<=25 and |b'|<=24*8+40=232, so |H_b|<=232.
The actual H=4u/(1+u^2) has |H|<=2 globally.

For each new sector choose the ordinary zero-mean Minkowski state
in the reference past and evolve from -3 to the anchor 0 with b.
The metrics and their operators agree on the open slab u>-1.
Restrict the prepared state to a smaller slab about zero, then
extend with the ACTUAL a-equations to all times. Uniqueness of
the linear Cauchy problem identifies this with the explicit formula

    C_added(0,k)=U_b(0,-3,k) C_flat(k) U_b(0,-3,k)^T,
    W_added(t,s,k)=U_a(t,0,k) C_added(0,k) U_a(s,0,k)^T.

In particular the actual past uses U_a, not U_b. This is neither
an alteration of the physical bounce nor a finite-order
instantaneous vacuum choice at its center.

For Proca, dE/du is the explicit derivative of the positive
Hamiltonian on solutions. Each term has scale power -3,-1 or 1,
so |E'|<=3|H|E. For the tensor use the positive modified matrix

    G_T=diag(a(1+q)/2,2/a^3).

Its evolution cross matrix is [[0,a^-2],[a^-2,0]].
The equality

    (a T^2/2+2 PT^2/a^3)/a-2 T PT/a^2
      =(T-2 PT/a^2)^2/2

and its reflected-sign version give |E_T'|<=(3|H|+1/a)E_T.
These bounds are uniform in k and apply in either time direction.

From -3 to 0, the Proca energy exponent is at most
3*232+3*2=702. The tensor adds at most three, giving 705.
Over this preparation a lies in [1,25]. Proca energy eigenvalues
lie between 1 and 10^6(1+q); tensor modified energy eigenvalues
lie between 2/25^3 and 25(1+q)/2. Therefore

    ||U_P,b(0,-3)|| <=1000 sqrt(1+q) exp(351),
    ||U_T,b(0,-3)|| <=(625/2) sqrt(1+q) exp(705/2).

These constants are intentionally loose but explicit, finite
and valid at every frequency, not just asymptotic modes.
