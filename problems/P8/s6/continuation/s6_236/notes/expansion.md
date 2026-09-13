# Complete diagram cancellation before estimation

Keep H0=ln(n)-Log(L-i0), with L=1-s xi(1-xi). The exact primitive jets and its full mass derivative are

T=(H0-1)/n+[3L H0+(b-7L)/2]/n²+R_T,
-partial_n T=(H0-2)/n²+[6L H0-10L+b]/n³+R_D.

Differentiating ln(n) supplies the minus1 in the first box term and the additional minus3L in the second. Treating this logarithm as a constant gives a wrong result.

Let B(s)=-integral Log(L-i0), J(s)=integral L Log(L-i0) and I(s)=integral L=1-s/6. After all parameter averages, the complete jets are

Cbar(s)=[ln(n)-1+B(s)]/n
       +[3I(s)ln(n)-3J(s)-(7/2)I(s)]/n²+R_C,

Dbar(s,t)=[ln(n)-2+B(s)]/n²
         +[6I(s)ln(n)-6J(s)-10I(s)+t/6]/n³+R_Dbar.

The term t/6 comes from the second, HEAVY channel. The first variable is not interchanged with it.

For A_s=C+g²/(n-s), the literal S233 contact gives the exact identity

A_s/g²=-2/n+(s-4)/n²+r_a,
r_a=n^-3[s²/(1-s/n)+4(4/n-1)/(1-2/n)²].

Insert the full jets into the COMPLETE S235 expression

A1_base=g4/(16pi²) sum_channels[
 (A_s/g²)² B(s)/2+2(A_s/g²) Cbar(s)+Dbar(s,t)+Dbar(s,u)].

For one channel, use t+u=4-s only after keeping every term. The entire n^-2 term is -2ln(n); the entire n^-3 term is

2(s-4)ln(n)+8/3-7s/6.

Every B(s) and J(s) coefficient is zero in these orders. This cancels their imaginary parts as well as their real logarithms. Omitting even one bubble or a box term invalidates the cancellation.

Now sum all three channels using s+t+u=4. The result is

A1_base=g4/(16pi²)[-6ln(n)/n²+(-16ln(n)+10/3)/n³]+R_base.

Both displayed terms are independent of the on-shell kinematics. They are NOT set to zero in the base scheme: S235 showed that its finely cancelled symmetric value needs a finite adjustment.

The ALREADY specified V2S-T1-OS4 contact is exactly -A1_base at s=t=u=4/3. Consequently the two universal terms cancel in the difference between a physical configuration and that symmetric point. No new value or derivative counterterm is introduced, and no frozen prescription is edited. The complete remainder proof, rather than power counting, controls everything left over.
