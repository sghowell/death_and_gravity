# Nonlinear metric map and the already fixed profile

Let vhat be the hat-metric log-scale source and retain the physical lapse. The actual physical log-scale is

    v = vhat + omega(u,N),
    omega = -1/4 log[(h-1+N^-2)/h],  h=(1+u²)^3.
    omega_N|1 = 1/(2h),
    omega_NN|1 = -3/(2h)+1/h².

Every needed time jet is differentiated with its background coefficient retained. Pulling back just the quadratic density along the linear map misses the first variation acting on the second map jet. After weighted integration by parts that extra quadratic contact is 3 p_local omega_NN n²/2.

An independent check composes the exact second Taylor coefficient of the full covariant density: its quadratic term receives the two first-map jets, while its linear term receives D_u^j(omega_NN n²/2). Both Euler outputs of the difference from the compact contact representation vanish at all three adiabatic orders. This avoids assuming a stationary metric action when the metric piece alone has a nonzero one-point function.

The S6.82 profile has already been chosen and frozen. By linear decomposition, its local matched contribution in the flat cutoff plateau is

    P_local = -p_local + (rho_local+p_local)(1-N^-2)/2.

Its actual quadratic density is obtained by differentiating N exp(3 omega+3 vhat) P_local with respect to N,vhat, then evaluating N=1,vhat=0. No state or counterterm is reselected under variation. Adding this contact to the pulled-back metric density cancels the constant m^4 term identically: P_local=-5/2 for that order and the metric density is +5/2.

At fourth derivative order the resulting two-source coordinate Euler coefficient is exactly

    -4 [[1/(4h²), 1/(2h)], [1/(2h), 1]]
    = -4 e e^T,  e=(1/(2h),1).

It is nonzero, negative semidefinite and rank one on I, with kernel (1,-1/(2h)). The pointwise fixed profile adds no derivatives and cannot change this top coefficient. No old invertible two-source prepared block, its inverse or its frozen pole transfers to this new rank-deficient local response.

