# Full generic triangle identity

Use R1+R2+R3+k=0, k^2=0 and epsilon*k=0, tr(epsilon)=0.
Keep independent symbols

    v=R1^2, w=R2^2, d=R1.R2, a=k.R1, b=k.R2,
    H=epsilon(R1,R1), J=epsilon(R1,R2), L=epsilon(R2,R2).

Then R3=-R1-R2-k, R3^2=v+w+2d+2a+2b,
k.R3=-a-b and epsilon(R3,R3)=H+2J+L. No scalar mass
shell or special physical Gram relation has been used.

Set x+y+z=1, M=x+y+n*z and U=x*y*R1^2+y*z*R2^2+x*z*R3^2.
For each of the three split lines, cyclically rotate both the vertices
and their original masses/weights. In the rotated frame
Q0=0,Q1=-R_first,Q2=-R_first-R_second. The split endpoints give
Q0_eff=gamma*k, with the other two Qs unchanged. Set
barQ=sum x_i Q_i and
U_gamma=sum_(i<j) x_i*x_j*(Q_i_eff-Q_j_eff)^2.

The exact S342 TT density is
-2*x_split*epsilon(Q0-barQ,Q0-barQ)/(M-lambda^2 U_gamma)^2.
The numerator carries lambda^2 under common momentum scaling.
Its degree2r coefficient is therefore

    -2*r*x_split*epsilon(Q0-barQ,Q0-barQ)
       integral_0^1 U_gamma^(r-1) dgamma / M^(r+1).

The generic calculation retains every line, every split weight and the
whole gamma polynomial. Subtracting the explicit labeled-box contact
at the SAME order gives zero at r=1,2 and, identically in every listed
Gram/TT variable and x,y, at r=3 gives

    -8*x^2*y^2*z^2*(a^2*L+b^2*H-2*a*b*J)/M^4.

This is a whole symbolic simplex identity, not an inference from
sampled kinematics or a Ward-only guess. Reversing the insertion sign
fails already at the lower order.

The bracket equals the literal linear-curvature contraction
R_abcd R1^a R2^b R1^c R2^d in the S336 convention. Replacing the pair
(R1,R2) by(R2,R3) changes its transverse edge vector only by-b*k,
so leaves the contraction unchanged. Thus the two singleton scalar
momenta may be used for the subsequent full Bose normalization.

Finally x=(1-z)*xi,y=(1-z)*(1-xi), with simplex Jacobian1-z.
The xi integral is integral xi^2(1-xi)^2=1/30, giving
-4*z^2*(1-z)^5/15. All three masses are retained in M; the exact
pointwise identity itself is valid for any positive mass assignment.

