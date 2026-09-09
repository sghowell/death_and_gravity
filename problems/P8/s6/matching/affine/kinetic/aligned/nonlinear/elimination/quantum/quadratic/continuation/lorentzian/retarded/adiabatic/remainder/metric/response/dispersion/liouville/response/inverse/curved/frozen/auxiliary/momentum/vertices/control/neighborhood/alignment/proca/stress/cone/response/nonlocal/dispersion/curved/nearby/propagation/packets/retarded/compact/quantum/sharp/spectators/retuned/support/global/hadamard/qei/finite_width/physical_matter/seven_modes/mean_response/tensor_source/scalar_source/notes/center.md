# Independent center jets and a local bounce representative

Let each source be represented by its natural-phase Hessian G.
The exact state transport derivative is

    D G=G'+M_natural^T G+G M_natural.

At u=0 the selected addition is I4, so its source is
tr G/2. Applying D again gives its actual derivatives;
no frozen-frequency or minimally coupled equation is used.

For q=k^2 the anchor sources per unit added amplitude are

    F_N=(542720000q^2-672623400q-110222113)/393660000,
    F_xi=0,
    F_p=(10880000q^2+33455400q-6985738267)/196830000,
    F_psi=(40000q^2+487200q+5993299)/9841500.

The source derivatives and actual mean coefficients give

    n=F_N/(2J),
    dp'=ell beta n+F_p,
    xi''=-dp'/2+alpha' n/3+F_xi',
    n''=[2alpha' dp'-3ell beta xi''+F_N'']/(2J)-n J''/J.

Here J=243/160, J''=-633/200 and alpha'=9.
The proper Hubble derivative variation is xi''+n''/2-11n.
Every coefficient is retained in the report. Direct absolute
polynomial sums on 1<=k<=2 give

    |n|<10eta, |xi''|<40eta,
    |n''|<200eta, |delta Hdot_proper|<200eta.

The acceleration kernel has opposite signs at the two band
endpoints. This is not a claim about the sign of the chosen
band average. The magnitude bound does not depend on that sign.

On the strip define the exact-frame representative

    N_eta=1+n, hat_a_eta=a exp(xi),
    a_eta=hat_a_eta[(h-1+N_eta^-2)/h]^(-1/4).

For eta<=10^-20, |log(a_eta/a)|<=|xi|+|n|<130eta
and N_eta>1/2. Parity gives zero physical Hubble parameter
at the anchor. The exact second logarithmic derivative is

    4+xi''-(3/2)[(1+n)^2-1]+n''/[2(1+n)].

The absolute jet bounds put it above
4-300eta-150eta^2>3. Division by N_eta^2 leaves a
strict local bounce. This does not prove completeness of
this scalar-perturbed representative, global smallness,
a quantum state co-evolved on it, or an exact SEE solution.
