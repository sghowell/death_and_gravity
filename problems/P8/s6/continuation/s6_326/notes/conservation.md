# Conservation must be differentiated before it is imposed

The soft Lorentz vector field on leg i is
delta k_i=P_ii*q/D_i-eta*A*k_i. It preserves each k_i^2, but on the
conservation surface sum delta k_i=Sbar*q, not zero.
For a null leg j, Z_j=sum_(i!=j) k_i.k_j vanishes on that surface while
delta Z_j=D_j*Sbar. Equivalently sum_(i!=j) S_ij=D_j*Sbar.
Thus differentiating a kernel after deleting its conservation-zero terms
is unsafe.

The frozen S300 finite angular kernel K_f has massive pair contribution
f_signed(z)/2 and every pair involving a null leg contributes
2*z*ln(2*|z|/(w_i*w_j)), omitting the w factor on a unit massive leg.
Here f(t)=2*(2*t^2-1)*acosh(t)/sqrt(t^2-1), t>=1, and
f_signed(z)=sign(z)*f(|z|).
The auxiliary-mass logarithmic kernel is
K_raw=K_f-2*sum_null ln(m_j/w_j)*Z_j+o(1).
It has the same on-shell value in the limit but NOT the same derivative.
The real logarithmic coefficient is
-delta K_raw/(8*pi^2)+Sbar*sum_all D_j*ln(|D_j|/m_j)/(4*pi^2).
Its mass logarithms cancel using delta Z_j=D_j*Sbar. The finite answer
contains ln(d_j), d_j=D_j/w_j, not ln(D_j). Individual ln(w_j) cancel;
taking absolute values earlier would invent a multiplicity-sensitive
sum w_j*|ln(w_j)| loss.

## Finite coefficient

All following pair sums are UNORDERED. Write alpha_ij=|k_i.(1,n_j)|
for a massive/null pair, delta_ij=1-n_i.n_j for a null/null pair, and
c(t)=t*(2*t^2-3)/(t^2-1)^(3/2). Then C=F+i*G with

F = -sum_MM f'(|z_ij|)*S_ij/(16*pi^2)
    -sum_MN [ln(2*alpha_ij)+1]*S_ij/(4*pi^2)
    -sum_NN [ln(2*delta_ij)+1]*S_ij/(4*pi^2)
    +Sbar*[sum_M D_i*ln|D_i|+sum_N D_j*ln d_j]/(4*pi^2),

G = [sum_same_MM c(z_ij)*S_ij
      +2*sum_same_MN S_ij+2*sum_NN S_ij]/(8*pi)
    -Sbar*sum_outgoing D_i/(4*pi).

The same MN pairs have an outgoing massive leg. For complex A, F and G
are complex-linear components, not separately the real/imaginary parts
of the number C. The outgoing sum is exactly 2E. Detector-resolution
terms are not included in C and are not erased. The physical named term
is C*ln(1/omega)*H_tree/kappa^(3/2).
