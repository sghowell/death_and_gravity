# Cancellation-free on-shell anchors

Let A=x(1-x) and b_n=integral_0^1 A^n dx.
Integration gives b_n=(n!)^2/(2n+1)! and
b_(n+1)/b_n=(n+1)/[2(2n+3)].
Since 0<=A<=1/4, b_n<=4^-n.

Expand J(s)=sum_(n>=1) b_n s^n/(n mF^(2n)).
On the entire two-point domain this expansion and all
needed derivatives converge uniformly: |s|/(4mF^2)<1.
Grouping before evaluation gives, with C=2NY/Q,

    a_n=3b_n/[n(2n+3)] > 0,
    f(1)=f(0)+C[2/3-sum a_n/mF^(2n)],
    f'(1)=C[2/3-sum (n+1)a_n/mF^(2n)],
    f_R(0)=-C sum n a_n/mF^(2n).

The three coefficient regroupings are checked independently.
In particular f_R(0) is not recovered by floating-point
subtraction of f(0) and f(1), both of order 10^193 at
the selected boundary.

For n>=1, respectively bound a_n/b_n, (n+1)a_n/b_n
and n a_n/b_n by 3/5,6/5,3/5. Positive factorizations
in the source prove all three inequalities.
After retaining K terms put x=1/(4mF^2). The common
positive omitted sum is at most

    T=(3/5)x^(K+1)/(1-x).

Subtract T from the mass-increment lower endpoint and
2T from the slope lower endpoint; subtract T from the
negative zero-remainder lower endpoint. These enclose
the infinite series, not just the retained terms.
All reported endpoints use exact rational arithmetic.
The public interface accepts exact finite rational
mF>=36 and one to sixteen native-integer terms.
