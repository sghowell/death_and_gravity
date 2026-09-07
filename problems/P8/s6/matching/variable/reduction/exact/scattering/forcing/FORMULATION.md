# S6.36: fixed-pulse zero-data own-f response has separated original-Q readouts

Keep the unchanged action, original off-shell physical metric g=eta and
clock theta=T, and the exact positive isotropic own-f branch of S6.33.
The linear tensor equation and canonical normalization are those of the
recursively rebuilt S6.35 certificate. This is not the rolling coupled-
relative coefficient-80 problem or a full physical-g/clock solution.

## One fixed physical input and a specified causal state

Fix M,tau,eta>0 and write u=T/tau, L=1/100, a=L/2. For every
0<delta=c-2<=10^-6 solve the actual own-f tensor equation

    (k Q_u)_u + 2B1 b(Q-q)=0,
    k=b^3/N, Q(-L)=Q_u(-L)=0.

Choose any one delta-independent real C-infinity pulse q with
0<=q<=eta, compact support strictly inside (-3L/4,-5L/8), and
integral(q du)>=eta L/16. The same q and tau are used for every delta.
An explicit flat-step pulse with support [-95L/128,-81L/128], ramp
width L/64 and plateau width 5L/64 proves this class nonempty.

The prescribed q is a physical-metric tensor perturbation, not a separately
conserved matter stress or the free-chi field. Maintaining it may require
driving g. The theorem concerns the linear response, not an unproved
uniform nonlinear remainder at finite perturbation amplitude.

## Actual loading is proved, not prescribed as hidden data

For Y=sqrt(k)Q, the canonical equation is

    Y_uu + V_delta Y = sqrt(k) mass q,
    mass=2B1 N/b^2, V_delta=mass-(sqrt(k))_uu/sqrt(k).

The full pump and source weight are retained. The frozen whole-box
bounds imply, on the fixed negative interval [-L,-a],

    1<k<5, 0<V_delta<80044, sqrt(k) mass>1/L^2.

A positive-Green first-zero proof on source-to-load lags at most L/4
gives G_u>=59989/80000>2/3 and G/h>=219989/240000>9/10. Hence

    Y_delta(-a)>=219989eta/30720000>9eta/1280,
    Y_delta,u(-a)>=59989eta/12800>eta/(24L).

Only the relaxed area and support conditions are used in these bounds.
The stronger uniform constants remain before taking limits. On [-L,-a],
the denominator D>=2-2/(1+a^2)^4>0, so the actual coefficients and
zero-data response converge to punctured delta-zero limits. This does
not introduce a literal action at the singular center.

With r=sqrt(u^2+delta/8), t=asinh(sqrt(8)u/sqrt(delta)),
Y=sqrt(r)psi and Z=(psi,psi_t/rho), rho=sqrt(7)/2, the loaded state
w_delta=Z(-a) has one fixed real nonzero limit w. In particular
norm(w)>9eta/(1280sqrt(a)).

## Uniform finite-error propagation suffices

The pulse has ended before [-a,a]. Recomputing the S6.35 comparison at
this smaller half-width gives central integral<33/25000, both tails
<=3/800 and normalized budget<507/125000<1/200. Thus the actual
Z propagation differs in operator norm by less than epsilon=4/399
from the reference Cauchy map

    F_ref(theta)=V0 D(theta) T_t D(theta) V0^*,
    theta=rho asinh(sqrt(8)a/sqrt(delta)).

Both endpoint phases are retained. The reference's A_J part rotates
with 2theta while its B_J part is constant. Its oscillatory norm is
abs(A_J)>1. The error is uniform, not assumed to vanish.

For the scalar readout choose the fixed phase offset
theta0=arg[A_J(w1+i w2)]/2 modulo pi in [0,pi). The exact sequences

    delta_n^+=8a^2/sinh^2[(theta0+pi n)/rho],
    delta_n^-=8a^2/sinh^2[(theta0+pi n+pi/2)/rho], n>=2,

are positive, lie below 10^-6 and tend to zero. The first reference
components have separation 2abs(A_J)norm(w). Keeping BOTH finite
propagation errors leaves a separation greater than (790/399)norm(w).
The actual errors need not converge, and no actual subsequence is
asserted to tend to its reference center. The bound n>=2 concerns
domain membership, not a claimed finite-n separation rate.

## Readout in the original hidden metric, with its physical clock

The exact endpoint identity is Q(+a)=sqrt(r/k) Z_1(+a). Its prefactor
converges to sqrt(a/k0)>0, with k0<5<25/4. The full endpoint map
retains k_u and uses Q_u=tau Q_T; no metric or matter frame is changed.
The chosen phase pair therefore satisfies

    liminf_n [Q_delta_n^+(+a)-Q_delta_n^-(+a)]
      > (790/399)(9/1280)(2/5)eta
      =237eta/42560>eta/200.

The actual endpoint family is bounded, so in particular

    limsup_delta->0 Q_delta(+a)-liminf_delta->0 Q_delta(+a)>eta/200.

This is the unchanged original Q, not an adjustable projection of a
two-vector. Neither Q nor the actual (Q,Q_u) pair is Cauchy in the
delta limit. No separate derivative-readout oscillation bound is claimed.
The complete physical experiment runs from T=-tau/100 to T=tau/200;
its duration and its input pulse do not shrink with delta.

## Scope boundary

This excludes a delta-independent limiting endpoint response for the
stated fixed-input, zero-hidden-data experiment. It does not exclude
delta-dependent phase-retaining descriptions, other states or input
classes. Smooth compact pulses are not a proved low-frequency EFT band.
No cutoff, all-order reduction, conserved matter loading, full CD/matter
dictionary, physical vacuum, positivity/V/G verdict, UV completion or
original P8 closure follows. A causal inverse is not silently inserted
into a single-copy variational action.

The level is CERTIFIED by written continuous proofs, exact margins and
independent source/readout audits, not proof-assistant FORMALIZED.
