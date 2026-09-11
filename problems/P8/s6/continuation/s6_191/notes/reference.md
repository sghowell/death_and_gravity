# Parameter derivatives of the finite ordered reference

Keep the actual full S6.190 reference rhat=sum_(n=1..10)r_n,
not a new adiabatic state. Differentiate its ordered
Riccati recursion through amplitude order2. For mixed
index beta=(d,a), distribute derivatives among factors
with the complete product-binomial and multinomial
coefficients. The inverse-frequency mixed derivatives
are those in mixed.md. They already include derivatives
of omega itself.

The resulting explicit positive majorants satisfy

    ||partial_t^d partial_epsilon^a r_n||
       <=b_(n,d,a) omega^-n,  n+d<=11, a<=2.

No matrix product is commuted. In particular all
parameter derivatives preserve the reference's inverse
frequency order. The finite residual identity remains

    F=r10'-[R,r10]
      +sum_(j+l>=10, 1<=j,l<=10)r_j S r_l.

Differentiating this exact finite expression, not a
discarded formal tail, gives constants A_a,C_a with

    ||partial_epsilon^a rhat||<=A_a nu_minus^-1,
    ||partial_epsilon^a F||<=C_a nu_minus^-10.

The rational constants are fully recorded. Their
approximate sizes, only for orientation, are

    A0~5.732, A1~125.1, A2~6493;
    C0~3.623e23, C1~1.180e26, C2~4.233e28.

The code checks all zero-parameter coefficient jets
against the prior finite reference, not just one
sampled leading term. That includes the exact residual
and the mass-normalized reference norm. The same
frequency lower envelope and K=1e16 analysis partition
remain; no physical mode is discarded.
