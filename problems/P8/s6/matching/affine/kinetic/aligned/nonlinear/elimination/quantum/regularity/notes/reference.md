# Eighth-order reference without changing the state

Keep the S6.55 all-order Gaussian state, all its fixed cutoff functions,
the S6.56 action and the S6.53 finite prescription. Work in tau=1 units
on I=[-1/2,1/2], with m>=1000, a=(1+u^2)^2 in [1,A], A=25/16.
Write omega^2=m^2+k^2/a^2 and nu^2=m^2+k^2/A^2. Then
nu<=omega<=A*nu, and nu is independent of time at fixed comoving k.

The reference used to estimate the exact state is now

    W=omega*S, S=1+sum_{n=1}^4 P_n t^n, t=omega^-2.

The P_n are the actual, unchanged S6.55 all-order recurrence
coefficients. This reference is not a new state preparation.
Let D0=partial_u+z' partial_z, lambda=-H*z,
R=D S=sum B_n t^n, B_n=D0 P_n-2n*lambda*P_n. Since
D t=-2lambda*t, the exact residual f''+(omega^2-U)f=rho*f of
f=(2W)^(-1/2) exp(-i integral W) is

    rho=(1-S^2)/t+2P_1+T/S+3R^2/(4S^2),
    T=sum [-D0 B_n/2+(n+1/2)lambda B_n] t^n.

This follows by expanding W'/W=lambda+R/S and using
2P_1=-U-lambda'/2+lambda^2/4. Multiplying by S^2 gives a
finite polynomial. Its coefficients of t^0,...,t^3 vanish exactly
for both actual polarizations; the next term is O(t^4).

## Explicit differentiated residual bound

For every rational coefficient C(u,z), `box_bound` supplies both
|C| and |D0 C| on the full rectangle |u|<=1/2,0<=z<=1, with
an exact reconstruction. Polynomial multiplication propagates these
two envelopes by (V,D)*(W,E)=(VW,DW+VE). Thus cancellation of
the low coefficients is used exactly before bounding the remaining
ones, without guessing a high-order constant.

If (V_i,D_i) are the remaining numerator envelopes, put

    Vsum=sum_{i>=4} V_i m_min^[-2(i-4)],
    Dsum=sum_{i>=4} (D_i+4i V_i) m_min^[-2(i-4)],
    Rmax=sum |B_n| m_min^(-2n), m_min=1000.

The exact bounds verify |S-1|<1/2 and 2+2Rmax<4. Using
|lambda|<=2, |(S^-2)'|<=16|R| gives

    |rho|<=C/omega^8, C=ceil(4Vsum),
    |rho'|<=D/omega^8, D=ceil(4Dsum+16Vsum*Rmax).

Actual constants (transverse, longitudinal) are
C=(5207948696836,5165763544817) and
D=(473081396588352,474200988057393).
The report retains the rational coefficient envelopes as well.

The reference obeys omega/2<W<3omega/2 and |W'/W|<4.
For the physical p_f=f'-d f, |d|<=3, hence
|f|^2<=1/omega, |p_f|^2<3omega and |f p_f|<2.
For example 9/4+25/m_min^2<3 is a sufficient square bound.

## Oscillatory evolution with nonzero initial mixing

Write the exact prepared mode v=Acal*f+Bcal*conj(f), with the
same relation for first derivatives. Remove the diagonal phase
exactly as in S6.54: Acal=e^(-i eta)*a, Bcal=e^(i eta)*b,
eta'=r=rho/(2W), Phi'=Psi=W+r. Then

    a'=-i r e^(2i Phi) b, b'=i r e^(-2i Phi) a.

The initial coefficients here are not (1,0). The preparation
bound in notes/preparation.md gives |Acal_0|+|Bcal_0|<2.
The norm integral J<=2C/nu^9<1/4 gives exp(J)<2 and
|a|,|b|<4. The report checks
C/m_min^10<1/4 and (D+4C)/m_min^10<1, so

    Psi>=omega/4, |Psi'|<7omega,
    |g|<=2C/nu^10, |g'|<=(2D+64C)/nu^10,
    g=r/(2Psi).

Keep both endpoints and the initial coefficient:

    b(u)=b(u0)-[g*a*e^(-2i Phi)]_u0^u
                +integral (g*a)' e^(-2i Phi) du.

The endpoints contribute at most 16C/nu^10, the differentiated
g term at most (8D+256C)/nu^10 and the feedback term at most
8C^2/nu^19<2C/nu^10. Therefore, on the entire interval,

    |Bcal(u)|<=|Bcal_0|+Kmix/nu^10,
    Kmix=8D+300C.

A common Kmix is 5347035781757616. This estimate is on the
actual negative-frequency coefficient, not on an unremoved phase.
The exact CCR and the uniform preparation bound imply |Bcal|<1,
which justifies the quadratic-product estimates below.
