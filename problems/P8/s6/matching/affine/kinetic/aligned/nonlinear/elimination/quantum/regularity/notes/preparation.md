# Piecewise bounds for the frozen Cauchy preparation

S6.55 fixes W_B=W_4+sum_{n>=3} chi(nu/Lambda_n)P_n omega^(1-2n)
and its time derivative at u0=-1/2. The cutoffs are independent
of time, with Lambda_3=2m and

    Lambda_n=max(2Lambda_{n-1},2^n(1+C_n)),
    C_n=max(box|P_n|,box|D0 P_n+(1-2n)lambda P_n|).

Nothing in this continuation changes W_B, these cutoffs or the
exact subsequent oscillator solution. It only compares those
same data with the positive eighth-order reference W_8.

Let V_n,G_n denote the two displayed coefficient bounds. Besides
the previously computed orders three and four, the exact recurrence
and continuous boxes give

    V_5(T)=659109606040857/1024,
    G_5(T)=61205538081187995/1024,
    V_5(L)=640985109634975/1024,
    G_5(L)=60076898136335185/1024.

For an active coefficient with n>=4, omega>=nu>Lambda_n and
Lambda_n>=2^n(1+C_n). In any of the weights below, the
remaining denominator has power at least two. Each such tail
coefficient is bounded by 2^-n; the geometric sum is below one.
This applies to both values and time slopes because the cutoff
has no time derivative.

The three comparison bands admit
|W_B-W_8|<=F/omega^(r-1), |W_B'-W_8'|<=G/omega^(r-1):

| Band in nu | r | F | G |
|---|---|---|---|
| all nu>=m, used below 4m | 6 | 2V_3+1+V_4/m_min^2 | 2G_3+1+G_4/m_min^2 |
| nu>=4m | 8 | V_4+1 | G_4+1 |
| nu>=K | 10 | V_5+1 | G_5+1 |

Here K=max(8m,10^12). At m_min=1000, the exact frozen fourth
cutoffs obey 2Lambda_4<=10^12 for both modes. Since
Lambda_4=max(4m,16(1+C_4)), this same K completes their
fourth-coefficient cutoffs for every m>=1000. In the middle
band chi_3=1; in the high band chi_3=chi_4=1. No statement
about an unactivated infinite tail is needed.

For two positive Cauchy frequencies W,V, their half-log rates
J=W'/(2W), Q=V'/(2V) give exactly

    Bcal_0=(W-V+i(Q-J))/(2sqrt(WV)).

Use W=W_8>=omega/2, V=W_B>=omega/4 and |W'/W|<4.
Then |Q-J|<=2(G+4F)/omega^r, and consequently

    |Bcal_0|<=[2F+4(G+4F)/m_min]/nu^r.

Rounding up and taking both polarizations gives common constants

    B6=1813229, B8=906392614, B10=1536706170025.

The B6 envelope is valid for all momenta, and proves the initial
coefficient norm below two and the evolved |Bcal|<1. The
sharper B8 and B10 envelopes are used only in their justified
bands. The evolution contribution Kmix/nu^10 is integrated
separately over all momenta, preserving its stronger decay.
