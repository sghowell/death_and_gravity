# Constant Weyl²: physical tensor reduction and finite-window bounds

## 1. The new action and the fourth-order tensor equation

The candidate is S_CD/M1+cC integral sqrt(-g) C², with constant cC in the
original physical metric and beta=cC/M². The Bach-type first variation of
the extra term vanishes when the Weyl tensor vanishes; evaluate that full
variation before imposing FLRW. Hence the old background is exact. This
fact does not certify the candidate's scalar sector or interactions.

Conformal invariance of integral C² in four dimensions reduces its tensor
quadratic action to the flat-space expression. With TT polarization norm
E_ij E_ij=1 and a real Fourier component, the action is

    S0 = (M²/8) integral d eta a²(gamma'²-k² gamma²),
    S1 = (cC/2) integral d eta (D gamma)²,
    D = partial_eta²+k².

The independent action audit directly contracts the dynamic four-index
linearized curvature, performs the spatial Fourier average, identifies
the total derivative, and varies the literal cosmic-time action. Thus
the coefficient is not inferred merely from background C=0.

Variation gives

    E_eta = gamma''+2 calH gamma'+k² gamma
             - 4 beta D² gamma/a² = 0.

In cosmic time set q=k²/a², with q_dot=-2Hq, and define

    L0 = partial_t²+3H partial_t+q,
    B = Hddot-2H Hdot,
    J = 2(H²-Hdot),
    L1 = 8B partial_t-16Hdot q.

Direct operator multiplication, without imposing any field equation, gives

    a^-4 D² = (L0+J)L0 + 4Hdot q - 2B partial_t.

In particular, on L0 gamma=0,

    D² gamma = 4a² Hdot k² gamma - 2a³ B gamma'.

The exact fourth-order equation is Efull gamma=0, where

    Efull = L0 - 4 beta a^-4 D²
          = L0+beta L1-4 beta(L0+J)L0.

The chosen first-order reduced representative is Lred=L0+beta L1:

    gamma_ddot+(3H+8 beta B)gamma_dot+(1-16 beta Hdot)q gamma=0.

For its exact two-data solutions, the original-equation residual is exactly

    Efull[gamma_red] = 4 beta² (L0+J)L1 gamma_red.

This is a residual identity and formal first-order matching statement.
It does not bound a fourth-order inverse or establish a unique exact
four-data branch. Perturbative reduction must be distinguished from
retaining arbitrary extra higher-derivative solutions.
[Parker and Simon, sections I–II](https://arxiv.org/pdf/gr-qc/9211002).

## 2. A second route through the action, with the observable map

In cosmic time the Weyl term is (cC/2) integral a³(E0-2H gamma_dot)²,
where E0=L0 gamma. The off-shell first-order redefinition is

    gamma = y+2 beta E0[y]-8 beta H y_dot+O(beta²).

Its variation of S0 cancels both the E0² and E0*y_dot pieces of S1,
leaving the reduced action

    S_y = (M²/8) integral dt a³[(1+16 beta H²)y_dot²-q y²]
          + O(beta²).

The full map, including 2 beta E0, is needed for this off-shell equality
modulo a boundary term. On a perturbative branch E0[y]=O(beta), it reduces
to gamma=y-8 beta H y_dot+O(beta²). The first-order y equation is

    y_ddot+(3H+32 beta H Hdot)y_dot+(1-16 beta H²)q y=O(beta²).

Mapping back to gamma reproduces Lred. Using the y speed coefficient as
the physical metric's coefficient gives the wrong answer: at the bounce
H=0 but Hdot=4/tau², and the physical shift is -64 beta/tau² rather
than zero. Initial data and observables must also be mapped, not renamed.
The code tests the action cancellation and both first-order equation
coefficients independently. These are finite-order field-redefinition
identities, not a new physical matter frame.
[Solomon and Trodden, section II C](https://arxiv.org/pdf/1709.09695).

Solving a truncated reduced equation exactly is a choice of higher-order
representative. Other reductions can agree through first order yet differ
beyond it. Here that choice is explicit, and its original-equation residual
is retained; no all-order physical uniqueness is claimed.
[Glavan, introduction and section III E](https://arxiv.org/pdf/1710.01562).

## 3. CD coefficients, signs and derivative bounds

Use u=t/tau, ell=tau sqrt(1+u²) and x=u/sqrt(1+u²). Then

    ell² Hdot = 4-8x²,
    ell³ B = -56x+96x³,
    ell³ F = -448x+768x³,                     F=8B,
    delta c_T,red² = beta(128x²-64)/ell².

The sharp polynomial bounds are |ell²Hdot|<=4 and |ell³B|<=40. Thus
|delta c_T,red²|<=64|beta|/ell² and |beta F|<=320|beta|/ell³.
At |u|=1 the speed shift vanishes. Its sign reverses across those points,
so any nonzero constant beta produces a coefficient above one in one
region and below one in another. This concerns the specified physical
two-derivative representative only.

Cosmic differentiation of ell^(-w)f(x) acts as
ell^(-w-1)[(1-x²)f_x-wxf]. In particular

    ell^(n+1) H^(n) = 4(-1)^n n! T_(n+1)(x),

where T_j is the Chebyshev polynomial. Its absolute value is at most
4 n! for |x|<=1. The recurrence is checked through n=4, which is all
the residual calculation needs. The endpoints are compact limits and the
bounds hold at every finite centre time.

## 4. Fixed-unit old canonical norm and the retained k enhancement

Take the same physical gamma,gamma_dot data at t0, for both equations.
Set a(t0)=1, ell0=tau sqrt(1+u0²), and

    sigma=(t-t0)/ell0,     |sigma|<=1/100,
    Y=(M/2)a^(3/2)gamma,   P=dY/dsigma,
    w0=ell0² W0,          W0=q-(3/2)Hdot-(9/4)H²,
    N_t²=|P|²+w0|Y|²=2E0.

This P is ell0 times the physical-time momentum Y_dot. All norms below
use these fixed centre units. They are not unscaled physical-time energy
norms and do not choose a corrected quantum vacuum.

In physical time the same old Y obeys exactly

    Y_ddot+f Y_dot+(W0+DeltaW)Y=0,
    f=beta F,
    DeltaW=beta G,
    G=-16Hdot q-12HB.

The -12HB term is the volume-normalization contribution, not an optional
boundary convention. With q_dimless=ell²q,

    ell⁴ G = (128x²-64)q_dimless+2688x²-4608x⁴.

The non-q polynomial has maximum 392 and minimum -1920. Therefore
|G|<=64q/ell²+1920/ell⁴.

Require L<=k_com ell0<=4U with L=10^11 and U=10^12. From |ell_dot|<=1
and |H|<=4/ell, the window obeys ell/ell0 in [99/100,101/100] and
|log a|<=4/99. The looser a,ell/ell0 in [1/2,2] are sufficient. Hence

    L/2 <= ell0 sqrt(q) <= 8U,
    ell²q >= L²/16,
    W0 >= q/2 > 0.

The pinned tensor estimate |W0_dot|/W0<=18/ell follows also directly
from |ell³ mass_dot|<=132 and ell²q>=132. Over any subinterval of the
window its integrated bound is at most 4/11; the baseline propagator has
energy ratio at most 11/7 and energy-norm bound g=4/3. These internal
propagator bounds are used in the centre-data Duhamel argument.

Let epsilon0=cC/(M ell0)². The difference of the reduced and baseline
first-order generators, measured in the old energy norm, is bounded by

    ell0[|f|+|DeltaW|/sqrt(W0)]
      <= |epsilon0|[2560+4096U+122880/L]
      < 10^16 |epsilon0|.

For example 2560 comes from 320(ell0/ell)^3. The q term in DeltaW
becomes proportional to sqrt(q), not a frequency-independent constant,
after division by sqrt(W0). This is the explicit retained frequency
enhancement. Merely bounding |delta c_T²| would not bound an accumulated
phase difference.

## 5. Energy, Duhamel and first-Born remainder

Supply an independent constant bound
epsilon=|cC|/(M tau)²<=10^-30, so |epsilon0|<=epsilon. This also implies
|beta|q<=64U² epsilon<=6.4*10^-5 throughout the window. In particular
the physical speed coefficient and the y-action kinetic coefficient stay
positive under the stated, sufficient smallness condition.

Write C=10^16 and T=2/100 as a deliberately loose full-window duration.
For the reduced solution,

    |d_sigma E0| <= [18 ell0/ell
                     + 2 ell0|f|+ell0|DeltaW|/sqrt(W0)] E0.

The extra integrated exponent is at most C epsilon/25<=1/100.
The elementary enclosure exp(r)<=1/(1-r), for 0<=r<1, gives
exp(4/11+1/100)<=1100/689<5/3. Thus, relative to centre data,

    3/5 <= E0(t)/E0(t0) <= 5/3,

and the reduced solution's norm is at most g=4/3 times its initial norm.
The energy ratio excludes zero data; the following norm inequalities do
not. Duhamel and the old propagator bound give

    N_t(z_red-z_free) <= g² C T epsilon N_t0(z_initial)
                       < 10^15 epsilon N_t0(z_initial).

Let z_first be the derivative at epsilon0=0, with zero centre correction.
One further Duhamel insertion, using the same pointwise time-interval
estimate and integrating its length, gives

    N_t(z_red-z_free-epsilon0 z_first)
      <= (g³/2)(C T)² epsilon² N_t0(z_initial)
       < 10^30 epsilon² N_t0(z_initial).

The argument works in either time direction from the centre using absolute
interval lengths. It is the Taylor remainder of this second-order
representative, not an estimate of unknown higher-order EFT coefficients.
With |cC|<=1 independently and M tau=10^324, these become 10^-633 and
10^-1266. No claim of naturalness or optimized physical scale is made.

## 6. Bounded original-equation residual, without a singular inverse

In centre units, let h_j=ell0^(j+1)H^(j). The Chebyshev bounds and
ell/ell0>=1/2 give |h_j|<=2^(j+3)j!, j=0,...,4. Let qbar=ell0²q.
Evaluate (L0+J)L1 gamma on the exact reduced flow, retaining its extra
powers of epsilon0. The absolute coefficient majorants are

| Extra epsilon0 power | Coefficient of gamma | Coefficient of d_sigma gamma |
| --- | --- | --- |
| 0 | 102400 qbar | 1024000+10240 qbar |
| 1 | 27525120 qbar+65536 qbar² | 275251200+1310720 qbar |
| 2 | 1677721600 qbar | 16777216000 |

These follow by exact polynomial l1 arithmetic, not sampled times or
momenta. Because qbar>=1, |epsilon0|<=1 and |epsilon0|qbar<=1, divide
the gamma column by sqrt(qbar), weight the three rows by the extra
powers, and bound their sum by

    18760226816 qbar < 2*10^10 qbar.

For n_gamma²=|d_sigma gamma|²+qbar|gamma|², the old canonical change
gives (M/2)a^(3/2)n_gamma<=2N_t. The exact residual identity then gives

    ell0² |(M/2)a^(3/2) Efull[gamma_red]|
      <= 8 epsilon² (2*10^10) qbar N_t
       < 10^38 epsilon² N_t,

using qbar<=64U². This is dimensionally a bound in the fixed-unit norm
defined above. At |cC|<=1 and M tau=10^324 it is at most 10^-1258 N_t.
Nothing in this bound estimates the inverse of the singular fourth-order
operator or fixes its extra initial data. Those would require a distinct
physical-branch construction and error assumptions.

## 7. Causality and remaining matching obligations

In flat space the exact truncated fourth-order equation factorizes as
D(1-4 beta D)gamma=0. Its additional root lies at D=1/(4 beta), absent
from the reduced branch. Treating that root as a physical ghost or tachyon
of the EFT, or using the formal double metric-null principal symbol as a
UV front-velocity prediction, extrapolates beyond what is established.
The low-energy gamma coefficient and a formal infinite-frequency limit
are different questions. No full quantum superluminality verdict follows.

The cC term is a specified local candidate, not the full isolated matter
loop correction: cR, finite matching, anomaly, massless nonlocal and state
terms, other loops and higher operators can contribute. A scale-running
increment does not fix cC, and no spacetime-dependent scale is substituted
into this constant-coefficient action. Scalar/matter constraint reduction,
new cubic/quartic interactions, nonlinear stability, corrected quantum
initial states, backreaction and UV matching remain separate obligations.
The exact old background and a tiny conditional tensor response do not
certify those missing sectors or promote the old tree certificate to this
new action.
