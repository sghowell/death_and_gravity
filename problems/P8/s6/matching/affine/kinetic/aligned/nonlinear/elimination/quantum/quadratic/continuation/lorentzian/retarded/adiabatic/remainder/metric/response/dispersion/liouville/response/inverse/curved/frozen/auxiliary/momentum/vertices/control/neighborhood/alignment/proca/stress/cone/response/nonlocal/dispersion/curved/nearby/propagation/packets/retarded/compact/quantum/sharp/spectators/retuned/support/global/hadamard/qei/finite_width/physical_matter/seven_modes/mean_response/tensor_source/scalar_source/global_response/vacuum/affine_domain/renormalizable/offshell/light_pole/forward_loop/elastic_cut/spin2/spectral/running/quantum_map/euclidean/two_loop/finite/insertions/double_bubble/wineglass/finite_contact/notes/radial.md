# Two full heavy-mass integrals

For M>1 define

    J1=integral_0^infinity
             y/[(y+1)^2(y+M)] dy,
    J2=integral_0^infinity
             y/[(y+1)^2(y+M)^2] dy.

Their exact values are

    J1=M ln(M)/(M-1)^2-1/(M-1),
    J2=(M+1)ln(M)/(M-1)^3-2/(M-1)^2.

The code checks independent anchored primitives:

    P1=M[ln(y+1)-ln(y+M)]/(M-1)^2
           +1/[(M-1)(y+1)],

    P2=(M+1)[ln(y+1)-ln(y+M)]/(M-1)^3
           +M/[(M-1)^2(y+M)]
           +1/[(M-1)^2(y+1)].

Both vanish at infinity, and -Pj(0)=Jj.
The identity J2=-dJ1/dM is checked separately;
it is not the only justification for J2.

Consequently the same fixed contact is

    delta L_fin=6[C0 g J1+g^2 J2]/(16pi^2).

Both Jj are positive by their original
integrands, and J2<J1/M. For actual M>32,

    J1<2ln(M)/M<2ln(4M)/M.

This uses M^2/(M-1)^2<2 and discards the
negative second term in J1. The entire
integration domains extend to infinity.
No heavy inverse is expanded inside a loop.

The sign proof from the full Hessian gives
the sharper magnitude bound

    |delta L_fin|
       =6[|C0|gJ1-g^2J2]/(16pi^2)
       <6 LgJ1/(16pi^2),

where 0<|C0|<L. Estimating the two opposite
terms separately is unnecessary and weaker.
