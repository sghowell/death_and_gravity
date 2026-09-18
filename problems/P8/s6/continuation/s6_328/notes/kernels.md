# Jointly continuous single and double kernels

Write all outgoing conserved legs as massive k_i and positive null
w_j v_j, v_j=(1,n_j). Let D_i=k_i.q, d_j=v_j.q,
L_i=A(k_i,k_i)/D_i, ell_j=A(v_j,v_j)/d_j and
S_ab=D_b L_a+D_a L_b-2A(k_a,k_b), with null weights included
when referring to full legs. For any subset J,
sum_(a<b in J)S_ab=D_J L_J-A(K_J,K_J).
On the full conserved set this is zero. Thus the common ln2+1 in
the mixed and null pair kernels moves into a smooth massive term.

Define
F0=-sum_MM[f'(|k_i.k_l|)/4-(ln2+1)]S_il+S_M H_M,
where S_M=sum L_i and H_M=sum D_i ln|D_i|.
For a unit null v and alpha=|k_i.v|, set
K_M=L_i d ln(d/alpha)+D_i ell ln(|D_i|/alpha)
    +2A(k_i,v)ln alpha.
The physical massive gap alpha,|D_i|>=1/4 holds uniformly.
At d=0 define K_M=0. The unique joint limit is zero:
d ln d->0, alpha->|D_i|, |ell|<=2, and transversality gives
A(k_i,v)->0. No direction-dependent ell value is used.

For two unit null directions set delta=v.v' and
K_N=ell*d'ln d'+ell'*d ln d
    -[d'ell+d ell'-2A(v,v')]ln delta.
At delta=0,d>0 define K_N=2A(v,v)ln d. At d=0 or d'=0
define K_N=0, including their common intersection.
The TT Gram identity bounds the bracketed S by2delta.
As delta->0 away from d=0, S ln delta->0 and the first
two terms have the stated diagonal limit.
For d->0 with delta separated from0, use the equivalent expression
ell*d'ln(d'/delta)+ell'*d ln(d/delta)+2A(v,v')ln delta.
Here d'->delta, bounded ell multiplies a vanishing coefficient,
d ln d->0, and A(v,v')->0. The d' case is symmetric.
At the triple intersection use
|K_N|<=2|d ln d|+2|d'ln d'|+2delta|ln delta| ->0.
These cases cover all singular intersections, including joint paths.

The exact complete formula is
4pi^2 F=F0+sum_j w_j sum_i K_M(i;v_j,q)
             +(1/2)sum_j,l w_j w_l K_N(v_j,v_l;q).
The self diagonal j=l gives w_j^2 A(v_j,v_j)ln d_j from
the original Sbar_N H_N; it must not be omitted on the grounds
that there is no pair diagram with two identical labels.
All formulas are linear in complex TT A, and the proofs hold
uniformly over unit tensors and the compact massive physical domain.
At a collinear atom only kernel contributions vanish; its energy
and momentum still enter the massive recoil. Recomputing a smaller
radiative recoil state would be a different physical configuration.
