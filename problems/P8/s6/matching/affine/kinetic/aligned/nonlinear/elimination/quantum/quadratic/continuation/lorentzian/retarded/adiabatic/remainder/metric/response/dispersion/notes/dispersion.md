# Exact homogeneous massive dispersion

## Canonical current, contacts and subtractions

This is a flat, zero-temperature Gaussian comparator at fixed mass m>0.
Freeze the external clock's mass slopes alpha,beta and its second
slopes, but vary both physical N and logarithmic scale Z. The scalar
clock itself is fixed. This is not the rolling selected state, a
common-parent Minkowski vacuum, or an on-shell background assertion.

At the unperturbed point let omega^2=m^2+p^2 and z=p^2/omega^2.
For each physical oscillator J_X=-partial_X H has the form

    J_X=(A_X p_canonical^2+B_X omega^2 v^2)/2.

The independently differentiated Hamiltonian in S6.67 gives pair
vectors b_X=A_X-B_X, ordered (N,Z),

    b_T=(beta(1-z), 2(1-z)),
    b_L=(beta+alpha z, 2(1+z)).

With v=exp(-i omega t)/sqrt(2omega), Wick contraction gives
C_XY^+(t)=omega^2 b_X b_Y exp(-2i omega t)/8 per oscillator.
For R=i theta(C^+-C^-), its Laplace transform is

    omega^3 b b^T/[2(s^2+4omega^2)], Re(s)>0.

The full response also contains expectation(-H_,XY) at fixed physical
canonical variables. Differentiating that Hamiltonian twice, including
second mass slopes and N/scale factors, reproduces

    bubble(s=0)+contact=-partial_X partial_Y(omega(N,Z)/2).

This is checked independently for both sectors and all four entries.
Contacts are frequency independent and disappear from the thrice-
subtracted remainder, not from the full renormalized static response.

Define M(z)=2 b_T b_T^T+b_L b_L^T; the factor two is physical transverse
multiplicity. The exact fourth-order Taylor subtraction uses

    1/(4omega^2+s^2)-1/(4omega^2)
       +s^2/(4omega^2)^2-s^4/(4omega^2)^3
      =-s^6/[(4omega^2)^3(4omega^2+s^2)].

The 0,2,4 local coefficients belong to the existing prescription.
They are not set to zero by this identity. The two derivative orders
are independently matched to all 16 per-sector S6.67 adiabatic
coefficient entries; odd derivatives vanish in the flat stationary
comparator. Thus no new finite Wilson choice is introduced.

## Elementary radial integration

With measure p^2 dp/(2pi^2), y=p/omega and d=1+4m^2/s^2,
the subtracted bubble is

    R_sub(s)=-s^6/(256pi^2)
       integral_0^1 y^2 M(y^2)/[4m^2+s^2(1-y^2)] dy.

This integral is convergent at both ends for every real s>0.
At fixed s the original high-p tail is O(p^-3) in the radial
integral after the three subtractions. It is the massive threshold
formula, not an expansion valid only below threshold.

Write M(z)=M0+M1 z+M2 z^2. Then

    64pi^2 R_sub(s)/s^4=-(M0 I1(d)+M1 I2(d)+M2 I3(d))/4,
    In(d)=d^n atanh(1/sqrt(d))/sqrt(d)
           -sum_{k=0}^{n-1} d^k/[2(n-k)-1].

Polynomial division and a separately differentiated primitive check
every In. These real formulas are asserted for d>1. They do not
specify a complex continuation by independently choosing square-root
and atanh branches; a future contour calculation must continue the
original retarded integral consistently.

## Large positive frequency and the finite complement

For s/m tending to infinity,

    In(1+4m^2/s^2)
      =log(s/m)-sum_{k=1}^n 1/(2k-1)
          +O((m^2/s^2) log(s/m)).

This follows from atanh(1/sqrt(d))=asinh(s/(2m)); the remainder follows
by the ordinary smooth expansion of sqrt(1+4m^2/s^2) and log near one.
For n<=3 this is uniform over the compact clock-parameter interval.
It is an asymptotic statement, not a numerical error bound at s~m.

Let c=(alpha+beta)/4 and v=(c,1). Since M(1)=16 v v^T,

    64pi^2 R_sub(s)/s^4
       =-4 v v^T log(s/m)+C+O((m^2/s^2)log(s/m)),
    C=(M0+(4/3)M1+(23/15)M2)/4.

The leading radial pole is P=2 v v^T in the same 1/(64pi^2)
normalization. Both P and the finite local fourth coefficient F
independently replay S6.67 after alpha=4/(9h), beta=28/(81h).
The complete constant accompanying this logarithm is E=F+C, not F:

    F=[[1256/(6561h^2), -76/(243h)],[-76/(243h), -4]],
    E=[[38528/(98415h^2), 92/(135h)],[92/(135h), 14/15]].

The fourth-order isolated block is B4(s)=F+64pi^2 R_sub(s)/s^4.
This definition explicitly omits the tree and renormalized
zero/second-order terms. Tadpole contacts have no fourth derivative.

## Exact positive-real-axis inverse, with its limited meaning

Set Q=[[1,0],[-c,1]], so v^T Q=(0,1). Direct congruence gives

    Q^T M Q =
     [[16(1-z)^2/(2187h^2), -8(3-z)(1-z)/(81h)],
      [-8(3-z)(1-z)/(81h), 4(3z^2-2z+3)]].

Write Q^T B4 Q=[[A,B],[B,Cc]]. The integral above and d>1 give

    15616/(98415h^2) <= A <= 3128/(19683h^2),
    116/(243h) <= B <= 604/(1215h),
    Cc <= -4.

The first two endpoint limits are independently evaluated polynomial
integrals after cancellation of (1-y^2). For the last entry use
3z^2-2z+3=3(z-1/3)^2+8/3>0. The signs of the integrands establish
monotonicity between the listed endpoints. Therefore

    det(Q^T B4 Q)=A Cc-B^2 <= -62464/(98415h^2)<0.

Put a*=15616/(98415h^2), b*=604/(1215h). The absolute entries of the
inverse are bounded by [[1/a*,b*/(4a*)],[b*/(4a*),1/4]].
For 1<=h<=(5/4)^3, the maximum row sum is at most
1635582375/63963136<26. This is a matrix norm on the real axis,
in the Q chart, not a norm of a time-domain convolution.

As a control, det(E-2P log(s/m)) vanishes at
log(s/m)=-9137/58560. That is s<m, outside the large-frequency
regime. The strict exact determinant bound applies there too.
Consequently this asymptotic zero is not an exact positive-real
pole of the isolated block. Nothing here excludes complex poles
or establishes the full tree-plus-loop propagator's spectrum.
