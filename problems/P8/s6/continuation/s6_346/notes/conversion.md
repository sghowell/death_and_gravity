# Entire labeled-word conversion

Let W_(a,b,c)=(1/4)sum24 v_ij^a(p_k^2)^b(p_l^2)^c, a+b+c=3.
The labeled lift applies L=-Box_g powers to Phi^2,Phi,Phi.
For an endpoint momentum P with a_soft=P.k and H=epsilon(P,P),
its TT response is-2H*DD(v^r), where

    DD(v^r)=sum_(ell=1..r) binomial(r,ell)v^(r-ell)(2a_soft)^(ell-1).

This polynomial includes the continuous zero-transfer limit.
For each word, project its ENTIRE flat polynomial onto the rank6
S345 basis, then subtract that lift's full metric/connection response.
Use only sum p_i+k=0,k^2=0 and TT, with otherwise independent G,a,H.
The quotient by T is accepted only after the WHOLE residual vanishes.

    powers 003 012 021 030 102 111 120 201 210 300
    tau     -6  -2  -2  -6  -4   0  -4  -4  -4  -8

All ten flat projections and generic radiative identities are computed,
not inferred from one polarization. The complete three-channel v^3
bubble is half W_300, so its conversion coefficient is-4.

For the triangle, write U=xy*v+yz*w+xz*q, x+y+z=1, and expand
v^j U^(3-j). A multinomial (a,b,c) contributes light powers
xp=a+c,yp=a+b and heavy power zp=b+c. Its exact light integral is
xp!yp!/(xp+yp+1)! times(1-z)^(xp+yp+1)z^zp.
Multiply the coefficient tau_(j+a,b,c) and sum every term. This gives

    K0=-2(z-1)^4(55z^3+10z^2+4z+1)/35,
    K1=4(z-1)^3(6z^2+3z+1)/15,
    K2=-4(z-1)^2(2z+1)/3, K3=8(z-1).

The full covariance of the outer A(L) is essential: DD(AF)=A(v)DDF+
F(u)DDA, not F(v)DDA. S344 proves the lower curvature differences vanish,
so only A(0) multiplies its primitive degree6 difference.
