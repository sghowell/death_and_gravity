# S6.34: finite-parameter retarded own-f history

Keep the unchanged action and exact isotropic own-f background of S6.33,
including the original physical metric, clock and free canonical chi.
This is a linear tensor-response theorem on a prescribed off-shell
physical geometry. It does not impose the physical-g/clock equations.

The actual domain is M,tau>0, 0<delta=c-2<=1/625 and
x=T/(tau sqrt(delta)) in [-1/4,0]. The background parameters satisfy
v=delta x^2<=1/10000 and remain inside the frozen branch's whole box.
No literal c=2 action is introduced.

## Specified response and input class

With physical-metric anisotropy q prescribed and hidden anisotropy Q,
the exact finite-parameter own-f equation and state prescription are

    (k Q_x)_x+s(Q-q)=0,
    k=b^3/N, s=2Q_profile b/(D/delta),
    Q(-1/4)=Q_x(-1/4)=0.

Q_profile is the frozen beta1 numerator, not the hidden tensor amplitude.
The complete variable-kinetic derivative is retained. The physical clock
change and action measure are derived from the literal tensor action.

The test class contains the zero history and, for a fixed amplitude
eta>0, every C-infinity nonnegative pulse q of compact support inside
(-1/4,-1/8), with q<=eta and integral(q)>=eta/16. Positive response
lower bounds below refer to the pulse subclass. The zero history is
the comparison input with identical final metric germ and identical
zero initial hidden state.

## Finite-parameter coefficient and Green bounds

Exact inherited background enclosures and a separate polynomial proof give

    19/5<k<21/5, 42<s<65, 1<=D/delta<=3/2.

These are bounds for every admitted positive delta, not substitutes from
the limiting inner equation. A proof using positive rational coefficients
establishes the desingularized denominator identity throughout the domain.

Let R(x,a) be the retarded Sturm kernel with R(a,a)=0 and flux jump
k(a)R_x(a,a)=1. A first-zero argument proves positive flux throughout
0<=x-a<=1/4, followed by the quantitative bounds

    R(x,a)>=(7495/38304)(x-a)>=(3/16)(x-a),
    R_x(x,a)>=1415/12768>1/10.

The strict kernel comparison applies for positive duration; equality
holds at zero duration. A separate direct Volterra estimate gives the
weaker but sufficient ratio 46385/242592>3/16. A longer interval is not
silently certified when this sufficient positive-flux condition fails.

Variation of constants now proves, for every pulse in the class,

    Q(0)>=63eta/1024, Q_x(0)>=21eta/80.

The same proof bounds the entire window by
0<=Q<=975eta/2432 and 0<=Q_x<=325eta/152. An explicit smooth step/plateau
pulse proves the class nonempty without approximate integration; all
endpoint derivatives are controlled by an exact polynomial recurrence.

The pulse and zero input have identical entire final germs, not merely
matching finite Taylor coefficients. Their different hidden responses
exclude exact algebraic locking Q=q and any exact germ-local response
map for this stated test class. An arbitrary common output for the two
germs incurs Q(0) error at least 63eta/2048 on one input. Stateful and
nonlocal maps, or other choices of homogeneous data, are not excluded.

## Original physical-metric equation readout

Define the normalized physical-g equation by E_g=-2 delta S_T/delta q,
so E_g=M^2 q_TT+2beta1 b(q-Q). This is the same norm-two equation
convention as S6.33; the literal action derivative is -E_g/2. Its inner
normalization is

    Ehat_g=(tau^2 delta/M^2)E_g=q_xx+s(q-Q).

At the final zero-q germ, Ehat_g(0)=-(128/c)Q(0). Since
c<=1251/625, the pulse therefore satisfies

    abs(Ehat_g(0))>=(80000/1251)(63/1024)eta>3eta.

The physical equation magnitude includes M^2/(tau^2 delta), and the
literal variational derivative has half that magnitude. This is a
diagnostic of the original physical-metric equation, not a claim that
the prescribed history solves it without driving. No conserved matter
source or free-chi loading is assigned to this residual.

## Scope that the theorem does not settle

The physical window has duration tau sqrt(delta)/4. Its pulses vary on
a shrinking physical time scale; they are not a fixed physical low-
frequency class below an established heavy gap or EFT cutoff. The
amplitude bounds concern a linear response, not a nonlinear remainder.

The retarded inverse is defined by its differential equation and zero
data. Inserting a causal kernel into a single-copy action is not assumed
to produce a variational effective action. No all-state/inverse theorem,
full-parent solution, complete CD/M1 dictionary, scalar health, physical
cutoff, ultraviolet positivity/V/G result or original P8 closure follows.

The level is CERTIFIED by written continuous proofs, exact margins and
three separately authored scientific audits, not FORMALIZED.
