# Whole pure-gravity elastic finite cut, including evanescent terms

Let beta=sqrt(1-4mu/s), Q=s-4mu, and
Lg=log(Q/(4pi nu^2))+Gamma_E. A normalized massive bubble's phase is
pi*beta[1+eps(2+log(nu^2/Q))]+O(eps^2).
Combine eps=-EP, the raw loop factor and the exact normalized
sphere endpoint <1/(1-x)>=1/(2EP)+1.
This derives raw ImC0mumu=-pi/(s*beta)[1/EP+Lg].
The same massive cut of I4(t,s) gives
2pi/(s*beta*t)[-1/EP+log(4pi nu^2/(-t))-Gamma_E].

Write the original elastic tree as A=P/(1-z^2)+a0+a2P2(z).
The common box and massive triangle coefficients retain
VD=V+2mu^2 EP/(1+EP), HD=H-2mu^2 EP/[s(1+EP)],
with H=kappa(a0+a2P2). Expanding these coefficients times their
IR-pole masters supplies the finite trace terms; setting them to their
four-dimensional values before integration would lose them.

The full common-basis cut, normalized by beta/(32pi), is checked equal to

    P*A/EP + Faux + P*A*Lg + P1*A + P*A1,

where P1=8mu^2/(kappa Q), A1 is the complete S281 evanescent tree, and

    Faux=P^2 Fz-3P a2P2+a0^2+a2^2P2/5,
    Fz=[log((1-z)/2)/(1-z)+log((1+z)/2)/(1+z)]/2.

This is the entire pure-gravity specialization of S281, not merely its
leading soft logarithm. All finite angular terms are accounted for.
The massive bubble coefficient supplies precisely the remaining
polynomial angular piece after the shared box and triangle contributions.

As an independent algebraic calibration, summing all three massless and
massive bubble coefficients and then taking mu=0 gives
203(s^2+t^2+u^2)/40. This checks the whole crossed coefficient sum; it
does not license changing the physical external mass or asserting a
renormalized massive Wilson coefficient without rational matching.
