# Physical-source spectral cost

This child studies preparation of the frozen S6.23 regular tensor sector
in the unchanged S6.20 variable-coupling parent. It is not an additional
vacuum, source-free projection rule, local EFT or UV verdict. The companion
constructions and independently authored audits support the scoped
certificate below; original P8 remains open.

Use physical g proper time `u=T/tau`, `M,tau>0`, `c=2+delta`,
`0<delta<=10^-9`, `K=(tau*k_com)^2 in [1,4]`, and
`J=(a,b)=(-1/50,-1/100)`, of length ell=1/100. The physical source is a
real external conserved TT anisotropic stress Pi, with

```
sigma(u)=tau^2 Pi(tau*u)/M^2,
S_probe=(1/2) integral dT a_g^3 Pi gamma_g.
```

The actual four physical data `(g,g_u,f,f_u)` vanish at the left endpoint.
No homogeneous heavy data are reset after preparation. The Fourier norm
is that of sigma, not `j_H sigma` or a volume-weighted canonical source.
Its zero extension is used, with
`sigmahat(omega)=integral exp(-i omega u) sigma(u) du` and energy measure
`domega/(2pi)`. Physical temporal frequency is `omega/tau`.

On J the full physical operator has a regular delta=0 limit. This is a
punctured limiting operator, not a delta=0 physical action through the
bounce. Its exact physical g-source loading map is `C_K:L2(J)->R4`.
The output coordinates are `P0^-(b)^-1 Z(b)=Psi_-(|b|)^-1 Y(b)`, with
the actual canonical map and normalized outer variables of S6.21.
For `v in R2` the target is `Lv`, where `L=[I2;0]` is the regular-light
inclusion. Unit v is this stated normalization, not unit raw g Cauchy data.

For `0<=Omega<=100`, source budget S and target v, define the infimum

```
inf { E_>(sigma): sigma in H0^2(J), C_K sigma=Lv, ||sigma||2<=S },
E_>(sigma)=(1/(2pi)) integral_|omega|>Omega |sigmahat(omega)|^2 domega.
```

The bounded claims are:

1. Exact actual-source controllability, a zero-extended H0^2 Hermite right
   inverse, and a continuous uniform-K enclosure of the full Gramian.
   The continuously verified calibrated inputs are
   `9/(4*10^12) I <= C_K C_K* <=16/25 I`;
   `||sigma_H||2<2*10^8 ||x||`, `||sigma_H,uu||2<4*10^14 ||x||` for a
   nonzero full normalized target x. Every nonzero regular target requires
   `||sigma||2>19||v||`; every target with nonzero even component requires
   `||sigma||2>900|v_even|`. The zero target admits the zero source, with
   the corresponding non-strict zero bounds.
2. The exact L2 spectral optimizer and a one-parameter dual are characterized
   analytically. They are not numerically evaluated. For every budget
   strictly above the minimum L2 control norm, the H0^2 infimum equals
   the attained L2 minimum. Attainment in H0^2 is not automatic. A C5
   Hermite field construction gives a C1, H0^2 source, not a C-infinity
   source with strict interior support.
3. On the stated small band, the tail operator is at least `2I/3`.
   For unit regular v and S=10^6 the unconstrained spectral optimizer has
   `||sigma||2^2<=2*10^12/3<S^2`. Hence the budget is inactive, the exact
   relaxation has strict slack, and its infimum lies between the strictly
   positive lower bound `722/3` and the ceiling `4*10^12/9`. The even
   lower bound is `540000`. These are bounds, not a computed band optimum.
4. A sinc expansion through n=4 has continuous operator error at most
   `1/(3*11!)`, with the displayed loading-norm-squared multiplier for
   band Gramians. It reduces any later band-Gramian computation to nine
   actual loading moments. Those moments and a numerical band optimizer
   are not claimed to have been computed here.
5. A fixed source with `||sigma||2<=S` and limiting loading error e has
   actual final normalized error at most
   `42e+delta(8600||v||+12600000 S)` relative to Lv, or
   `42e+delta(8400||v||+12600000 S)` relative to the exact analytic prepared
   target. For unit v, S=10^6, e=0 and delta=10^-17, either is below 1/1000.
   This does not assert a small error at every delta<=10^-9 for that budget.

These are linear-probe statements, uniform on the declared spatial band.
The source may have either sign; its anisotropic-stress norm is not a
positive matter energy or a finite-amplitude backreaction bound. A Fourier
energy fraction is not a hard bandlimit, RMS derivative bound, rolling gap
or Wilsonian cutoff. In particular the entire small-band source energy
fraction is at most 1/3, so these sources are not claimed to be mostly in
that band. An unwindowed low-pass projection has nonlocal temporal tails;
it cannot be evolved with a newly imposed zero past state. Long preparation
outside the recorded local background remains a different problem.
No arbitrary-source local EFT, vacuum positivity, scalar/vector/nonlinear
health, original C/D operator matching or original P8 closure follows.
