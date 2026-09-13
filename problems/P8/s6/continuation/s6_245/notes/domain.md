# Whole scalar frequency and pair domain

All bounds refer to the full six scalar WKB iterates of S6.240. A four-jet extraction is used only to compare their UV coefficients, never to replace their physical values.

## Joint complex inverse-radius and time domain

Write \(P=|P|,\ U=m+P,\ \epsilon=1/100,\ x=1/r\). Allow \(|x|\le\epsilon/U\), \(|d-3|\le1/4\), and complex times within \(1/64\) of the real slab. At a real center \(a_0\le25/16\). The logarithmic derivative bound \(|2H|<10\) gives the relative inverse-scale-square allowance
\[
 e_a=\frac{10/64}{1-10/64}=\frac5{27}.
\]
The second-leg momentum allowance is \(e_g=2\epsilon+\epsilon^2\). Thus the full scaled squared frequency differs from \(a_0^{-2}\), relatively, by at most
\[
 e_\omega=e_a+(1+e_a)e_g+(25/16)^2\epsilon^2
 =\frac{14463467}{69120000}<\frac14.
\]
The branch continued from positive real frequency has
\(\Re(x\omega)>\sqrt{3/4}/(25/16)>11/20\) and \(|x\omega|<2\). Also \(|a|<100/59<2,\ |a^{-1}|<64/59<2\).

Set
\[
 O=a x\omega=\sqrt{1-2Pxu+P^2x^2+m^2a^2x^2}
\]
for the second leg, omitting the transfer terms for the first. Its full squared defect is bounded by
\(e_g+4\epsilon^2<1/40\). Consequently \(\Re O>49/50,\ |O|<11/10\), and the full mass fraction obeys \(|z-1|<1/100\). These statements include complex \(x\); the functions \(x\omega\) and \(O\) are the analytic normalized branches through \(x=0\).

## Every full iterate and readout

The scalar potential is \(U_d=dH'/2+d^2H^2/4\); \(|H|<5,\ |H'|<12,\ |d|\le13/4\) give \(|U_d|<100\). The complete base logarithmic derivative has \(|\lambda|<10,\ |\lambda'|<400\), and its residual bracket is below 400.

Apply the normalized S6.240 iteration to \(w=W/\omega\) on \(|w-1|\le1/100\). On each time loss \(1/1024\), \(|\log w|<1/50\). The first two Cauchy derivatives bound the entire nonlinear bracket by \(10^5\), including the squared logarithmic derivative; a \(10^7\) Lipschitz bound suffices. The coarser common allowance \(10^8\) therefore gives
\[
 q=\frac{4\cdot10^8\epsilon^2}{(10^{98})^2},\qquad
 |w_6-1|\le2q<10^{-3}.
\]
This is an induction on all six full analytic iterates, with no discarded higher-order terms. The frequency-square root stays in the same nondegenerate branch at each step. After six losses and one readout loss the remaining time width is
\[
 1/64-7/1024=9/1024>2R,\qquad R=1/20000.
\]

Let \(\bar W=a xW_6\). Then \(\Re\bar W>9/10,\ |\bar W|<2\). A Cauchy derivative gives \(|\bar W'|<2048\). The entire derivative feature
\[
 \bar d=\frac{(d-1)H+\bar W'/\bar W}{2},\qquad
 \bar\pi=-i\bar W-a x\bar d
\]
satisfies \(|\bar d|<5000,\ |\bar\pi|<3\). Moreover
\(\bar g=1/(xW_k+xW_l)\) has modulus below 1, since the sum has real part greater than 1.

## Complete normalized spatial pair

Let the five scalar Hamiltonian features be momentum, three gradients and mass times field. At \(N=1,Q=0\), the symmetric spatial vertex has blocks
\[
 M_D=\operatorname{diag}(-\operatorname{tr}D/2,\,
           (\operatorname{tr}D)I/2-D,\,
           \operatorname{tr}D/2).
\]
It obeys \(\|M_D\|_{\rm op}<3\|D\|_F\). Write \(\ell=-k+P e\), \(n=k/r\), and \(\bar\ell=-n+Px e\). With phase and \(a^{-3/2}\) removed, the complete normalized pair is
\[
 \beta_D=\frac{(\bar\pi_k,i n,m a x)^tM_D
                  (\bar\pi_\ell,i\bar\ell,m a x)}
                 {2a\sqrt{\bar W_k\bar W_\ell}}
       =\frac{a^3V_D}{r}.
\]
The gradient directions are not split into polarizations. The trace, gradient and mass pieces together have numerator at most \(14\|D\|_F\): use \(|\operatorname{tr}D|<2\|D\|_F\), \(|n\cdot\bar\ell|<2\), \(|n^tD\bar\ell|<2\|D\|_F\), \(|\bar\pi_k\bar\pi_\ell|<9\) and \(|m^2a^2x^2|<1\). The full normalization is below 2. Hence
\[
 |\beta_D|<28\|D\|_F<1000\|D\|_F.
\]
This bounds the whole tensor before the four angular geometric products are separated.

For all real internal momenta, the corresponding complex-time feature bound follows from \(W_6<9\Omega,\ 1/W_6<3/\Omega,\ |\partial_t\log(a^3W_6)|/2<20000\), and \(|k/a|<9\Omega\). The full canonical normalized feature vector is below \(25\sqrt\Omega\). Thus, on the required discs,
\[
 |a^3V_D|\le C\sqrt{\nu\mu}\|D\|_F,\quad
 |g|\le2/(\nu+\mu),\quad C=10^9,
\]
where \(\nu=\Omega_k,\mu=\Omega_\ell\). The scale factors cancel between the full pair and the volume; no nine-vector-polarization factor is present.

## Exact full-frequency UV identification

The direct recurrence for each normalized four-jet is
\[
 \bar W_{\rm new}^2=O^2+a^2x^2
 \left[-U_d-\tfrac12\partial_t{\rm rate}
                   +\tfrac14{\rm rate}^2\right],
 \qquad {\rm rate}=\bar W'/\bar W-H.
\]
The code computes all seven rows, initial plus six iterates, for both legs. The first iterate differs from the second. The second and sixth agree through degree four and all five coefficients equal the frozen S6.243 frequency. Thus its 140 complete endpoint coefficients are coefficients of the full family bounded here, not a fitted replacement.
