# Continuous logarithmic-tail estimates

For x>=1 set r=arcosh(1+x/2), v=tanh(r), and q=tanh(r/2).
Then sqrt(z)=v and

    exp(r)=1+x/2+sqrt(x(x+4))/2>=1+x,
    r>=log(1+x)>=log(2)>1/2.

The last strict inequality follows by integrating
1/t>1/2 on 1<=t<2. The same density becomes

    A(v)=16/15+v^2-3v^4,
    B(v)=v(3-2v^2+3v^4),
    D=A+B r, U=B/2, a=U/(D^2+pi^2 U^2).

The exact completions

    A+14/15=(1-v^2)(3v^2+2)>=0,
    B'=15(v^2-1/5)^2+12/5>0,
    B-(8/3)v=3v(v^2-1/3)^2>=0

give A>=-14/15 and 0<=B<=4. For r>=1, tanh(r)>3/4:
the positive exponential series gives
exp(2)>sum_(j=0)^5 2^j/j!=109/15>7.
Hence D>=2r-14/15>=r there. For 1/2<=r<=1,
the common threshold gap D>=16/15 also implies D>=r.
The same lower bound thus holds across the entire tail.

On 0<=v<=1 the literal polynomial derivatives give

    |A_v|<=14, |A_vv|<=38,
    |B_v|<=24, |B_vv|<=72.

Since |v_r|<=1 and |v_rr|<=2,

    |A_r|<=14, |A_rr|<=66,
    |B_r|<=24, |B_rr|<=120,
    |D_r|<=18+24r<=60r,
    |D_rr|<=114+120r<=348r,
    |U|<=2, |U_r|<=12, |U_rr|<=60.

Let Q=D^2+pi^2 U^2. Using D>=r, r>=1/2 and pi^2<10,

    |Q_r|/Q <=2*60+2*10*2*12*4=2040,
    |Q_rr|/Q <=2*60^2+2*348
                +2*10*(12^2+2*60)*4=29016.

The quotient derivative formulas now yield

    |a|<=2/r^2,
    |a_r|<=(12+2*2040)/r^2=4092/r^2,
    |a_rr|<=(60+24*2040+2*29016+4*2040^2)/r^2
            =16753452/r^2.

The x/r change of variables is essential:

    x a_x=q a_r,
    x^2 a_xx=q^2 a_rr+q(q_r-1)a_r,
    0<=q<=1, q_r=(1-q^2)/2.

Thus |q(q_r-1)|<=1, and the second scaled derivative is
bounded by 16757544/r^2. With C_tail=20000000 we have,
uniformly for every x>=1,

    max(|a|,|x a'|,|x^2 a''|)
       <=C_tail/r^2<=C_tail/[log(1+x)]^2.

This is a continuous differentiable tail bound, not a bare
asymptotic statement or a finite frequency truncation.
