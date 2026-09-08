# Actual-clock local envelopes and finite energy/pressure estimate

Use the [explicit prescription](prescription.md), the unchanged
S6.50 prepared modes and the S6.51 subtraction on I=[-1/2,1/2].
Set mu=m0, m=m0*tau>=1000 and first use tau=1 units. Relative to
1/(64pi^2), the local finite energy and pressure have the form
m^4*E0+m^2*E2+E4 and m^4*P0+m^2*P2+P4.

Direct substitution of the actual clock dictionary gives

    E0=-5/2-22/[27(1+u^2)^3], P0=5/2,
    E2=-32(1215u^8+3645u^6+3645u^4+1136u^2+7)
                   /[243(1+u^2)^5],
    P2=80(5u^2+1)/[3(1+u^2)^2],
    E4=-32(1701u^10+3969u^8+1782u^6-1009u^4-869u^2+22)
                   /[81(1+u^2)^7],
    P4=64(7u^4-1)/(1+u^2)^4.

These are first-variation coefficients in this prescription, not
estimates of unknown finite operators in a different matching model.
The original lapse jets and physical Hubble/time derivatives are
substituted before the bound. E0 has |E0|<=179/54<4 since h>=1.
For the others apply the S6.50 exact polynomial-box envelope to
P(u)/(1+u^2)^n: sum_i abs(c_i)/2^i bounds it on the entire I.
Exact polynomial reconstruction is checked. The resulting bounds are

| coefficient | energy absolute envelope | pressure absolute envelope |
|---|---:|---:|
| order zero | 4 | 5/2 |
| order two | 49537/648 | 60 |
| order four | 355657/2592 | 92 |

No time grid or unspecified order-one coefficient is used. Define
U_E,U_P by these respective triples. Since pi>3, their physical
local contributions relative to M^2/tau^2 are at most

    [U0*R^4+U2*R^2+U4]/[576*L^2],
    L=M*tau, R=m0*tau.

The [complex-dimension limiting proof](dimension-limit.md) identifies
the finite subtracted integral in the regulated prescription with the
physical S6.51 integral. Its uniform all-momentum bound is essential:
one cannot simply add a four-dimensional result to a dimensionally
continued local coefficient without justifying this limit.

The physical finite subtracted mode integrals are bounded by S6.51.
Adding the applicable local bound to each gives an explicit bound on
the vector's matched first-variation value in this fixed prescription.
At L=10^12,R=1000 the total energy and pressure bounds are

    51841169900859757/7464960000000000000000000000000,
    3125106104119/720000000000000000000000000,

respectively; each is less than 10^-14. The exact fractions, domain
and error decomposition are retained in the report. No sign of the
quantum correction is assumed in this absolute estimate.

This bounds the stated one-loop Gaussian contribution at fixed light
background. It is not a bound on other fields, unknown UV matching
coefficients, derivatives of the state term, the full scalar equation,
higher loops or the corrected constrained principal symbol. A small
unsigned correction does not protect a classically saturated matter
cone. State admissibility and a corrected solution remain separate
research; original P8 is not closed.
