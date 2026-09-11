# Actual chart and nonstationary retarded mean currents

The current exact map is

    g_physical=R^(-1/2)g_hat+(1-R^(-1/2))du du/X.

Here R^(-1/2) means the reciprocal square root, NOT R minus1/2.
The physical matter metric is unchanged. In u=t gauge the lapse
is unchanged and v=vhat+omega(N,t), omega=-log R(t,N^-2)/4.
The literal complete analytic R has exact clock jets
R=1, R_X=1/h, R_XX=0, h=(1+t^2)^3. Consequently

    w1=omega_N=1/(2h),
    w2=omega_NN=-3/(2h)+1/h^2.

This is differentiation at the reference clock, not replacement
of the full function by its clock polynomial at finite amplitude.

Let rho,P and delta rho,delta P be the full conditional physical
reference stress and its retarded variation. Define the two mean
Euler currents per FIXED reference volume a0^3. On the background
the physical currents are e=(-rho,3P), and their first variations
for physical sources (n,v) are

    delta e_N=-delta rho-3rho v,
    delta e_v=3delta P+3P(n+3v).

The volume/readout contacts are essential. The nonlinear map has
J=[[1,0],[w1,1]], with physical v=vhat+w1 n at first order.
Its second derivative adds the uncancelled term (3P w2 n,0).
Thus

    delta e_hat = J^T diag(-1,3) (delta rho,delta P)^T
                  + C (n,vhat)^T,

    C11=3P(w1+3w1^2+w2)-3rho w1,
    C12=-3rho+9w1 P,
    C21=3P+9w1 P,
    C22=9P.

The background vector stress is not zero; its second-map contact
cannot be discarded. No old scalar profile is added to cancel it.
An independent generic nonstationary local density checks this
entire chain rule. For a constant density c, rho=-c,P=c and the
physical stress derivative is zero, but at the bounce its actual
clock-current matrix is c*[[15/4,15/2],[15/2,9]].
The second-map term alone is -3c/2 in the NN entry.

For a general prepared state the response is a retarded mean-current
kernel, not an ordinary symmetric functional Hessian. Pointwise
input/output multiplication and the displayed local contacts
preserve causal support. The local-density Hessian is an independent
chain-rule control, not a replacement for the retarded kernel.
