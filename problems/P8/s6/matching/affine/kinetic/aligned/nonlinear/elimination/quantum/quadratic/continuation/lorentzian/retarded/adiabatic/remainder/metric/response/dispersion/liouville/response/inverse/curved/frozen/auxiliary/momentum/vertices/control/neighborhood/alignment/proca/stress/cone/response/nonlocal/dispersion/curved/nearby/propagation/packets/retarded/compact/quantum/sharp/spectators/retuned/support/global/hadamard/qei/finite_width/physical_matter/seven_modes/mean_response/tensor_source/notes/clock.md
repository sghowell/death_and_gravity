# Actual clock source and absolutely convergent transfer

Before clock restriction, X=phi_dot^2/N^2. Varying the
literal tensor Lagrangian with respect to phi_dot gives,
at phi_dot=N=1,

    partial L_T/partial phi_dot=-a^3 K/h,  K=r-3s.

The explicit h derivative vanishes only after this variation.
Consequently the Euler source in the convention delta S/delta phi
is E_phi=(K/h)'+3H K/h. The actual metric-source identity is
rho_T'+3H(rho_T+s)=-E_phi. Setting R=1 prematurely
would erase this source.

This is not merely a possible nonzero term: the selected
two-polarization band has E_phi(0)=0 and

    E_phi'(0)/eta=72+14q-4q^2

before averaging. Its lower bound 64 and upper bound 337/4
follow from the exact factorizations

    E_phi'/eta-64=2(2q+1)(4-q),
    337/4-E_phi'/eta=4(q-7/4)^2.

To control the entire integral put x=P^2/a^6,
y=qT^2/(4a^2), c=qTP/a^5 and omega=sqrt(q)/a.
The exact two squares r +/- c/omega prove |c|<=omega r,
and direct Hamilton evolution gives

    K'=4c+12Hx-4Hy,
    |K'|<=(4omega+12|H|)r,  |K|<=2r.

These pointwise quadratic inequalities persist under the
positive covariance and band average. Write t=1+u^2.
Since h'/h=3H/2, omega<=2/t^2 and r<=4eta/t^4,

    |a^3 E_phi|<=4eta [8/t^3+60|u|/t^2].

S6.104's exact integral bounds give the whole-line L1 upper
bound 1944eta/7<280eta. No conditional cancellation is needed
to define the clock exchange.

Finally a^3 E_phi=(a^3 K/h)' exactly and
|a^3 K/h|<=8eta/t tends to zero on both tails.
The signed whole-line transfer is therefore zero.
If mu=integral b(k)|k|^2 d^3k/(2pi)^3, then 1<mu<4 and

    integral_0^infinity a^3 E_phi du=(4-mu)eta.

For eta>0 the future transfer is strictly between 0 and 3eta;
the past transfer is its opposite. This is exchange with the
clock, not separately conserved tensor metric stress.
