# Exact physical pressure and its off-diagonal remainder

In the original helicity basis the spatially isotropic
pressure operator is p sigma1/3. The physical mass
rotation gives

    Pi=p^2 sigma3/(3omega)+pM sigma1/(3omega).

Its norm is p/3, not omega. For a negative occupied
projector with Bloch vector n, the energy and pressure
per helicity are -omega n_z and
-p^2 n_z/(3omega)-pM n_x/(3omega).
Neither is a rotating-frame frequency eigenvalue.

The four exact frames give

    n_x=cos(theta1)sin(theta3)
        +sin(theta1)cos(theta2)cos(theta3).

The first x-axis rotation does not change this component.
Put q0=pMdot/(2omega^3), b=q0dot/(2omega). The exact
next ratio is

    q1=b/(1+q0^2)^(3/2), theta1=atan(q1),
    b=p[Mddot-3M Mdot^2/omega^2]/(4omega^4).

The exact connection and first frequency correction
are retained. For the inherited bounds A=2pDelta/(tau E^3),
r=1024/(tau E)<=1/1024, we have |q0|<=A<1/4,
|q1|<=Ar and |theta_j|<=A r^j.

The required inequalities are

    |n_x-sin(theta1)|<=|theta3|
       +|theta1|(theta2^2+theta3^2)/2,
    |sin(theta1)-theta1|<=|theta1|^3/6,
    |atan(q1)-q1|<=|q1|^3/3,
    |b-q1|<=2|q1|q0^2.

The last bound follows from the derivative of
(1+u)^(3/2) on0<=u<=1/16. Summing gives

    |n_x-b|<=A r^3
       +A^3 r[2+(r^2+r^4+r^6)/2]
       <A r^3+3A^3r.

Use the S6.168 bound on n_z-(1-q0^2/2). Both helicities'
UV pressure terms are therefore

    P0=-2p^2/(3omega),
    P2=-p^2 M Mddot/(6omega^5)
       +p^2(p^2+6M^2)Mdot^2/(12omega^7).

The exact-state versus fourth-frame pressure error is
at most one third of its energy allowance because
||Pi||=p/3<=omega/3. The z-component frame and Taylor
errors obey the same one-third comparison. No
inequality is differentiated to infer pressure.
