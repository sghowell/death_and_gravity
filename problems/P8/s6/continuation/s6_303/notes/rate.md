# Uniform real known interference at all physical angles

Use the complete known finite S302 reference, not the matched
three-coordinate family. By Bose symmetry choose
tau=min(-t,-u), so0<tau<=q/2. Split the domain only for the proof,
not as an imposed physical angular cutoff.

## Near either endpoint,0<tau<=1

Relabelt=-tau,u=4-s+tau. Thenu<=-5/4 and all non-small channels
have magnitude at least1. For any negativex in[-12,0],
J_x<=1,J'_x<=1/6,|V_x|<=322 and|V'_x|<=28. Hence

|(V^2 J)'|<=105938/3<40000,
|(V J)'|<=245/3<82.

Forf=V^2J andg=VJ, the physical endpoint values
f(s)+f(4-s) andg(s)+g(4-s) are imaginary.
The REAL part of the two small-a boxes is therefore bounded by

40000(3+|ln tau|)+328.

This follows from the real spacelike divided differences; there is
no absolute1/t bound on the imaginary phase.
SimilarlyReB_soft(0)=0 and|d_t ReB_soft|<=82.
The complete finite soft division4(1/s+1/t+1/u)B_soft has
absolute real part at most800328, using|B_soft|<100000 for its
two regular denominators.

Use S302's full master and coefficient bounds for the other pieces.
The four other ordered boxes contribute at most
4*(110000*200*10+1300*200).
The two non-small nonbox rows contribute at most
2*(210000000*10+200000*1200+2500*100
   +2300000*12+30000000*21).
The three old finite Gram terms contribute less than22000.

The small-channel massless coefficients obey|c00|<=210000000 and
|b00|<=2300000. For0<tau<=1 their masters obey
|C_t|<=5/sqrt(tau)+|ln tau|/2+1,
|B00_t,finite|<=5+|ln tau|.

Add the entire regularized massive/rational bound from the previous
note. The coefficients of1,|ln tau|,tau^-1/2 in the resulting bound are

82760937423943327/1800,
107340000,
1050000000.

Each is below10^14. Thus

|ReF_known|<10^14[1+|ln tau|+tau^-1/2].

## Positive Born normalization

The unchanged S297 functions give
A_G>p/(2kappa*w),w=4tau(q-tau)/q^2,p=4q+16+8/q.
Sincew<=4tau/q andpq/8=V/2>=257/32>8,

A_m+A_G>=A_G>8/(kappa*tau).

Consequently the absolute known one-loop interference is at most

tau*|ReF_known|/(64pi^2*kappa).

For0<tau<=1,tau*|ln tau|<=1/e<1 and
tau+sqrt(tau)<=2. Withpi^2>9, this is strictly below
3*10^14/(576kappa)<10^12/kappa.

## Interior and complete coverage

If tau>=1, S302 applies directly withdelta=1 and gives
|A_known|<10^9/kappa^2. Also p/2>=257/18>14 andw<=1,
soA_G>14/kappa. The relative interference is below
2*10^9/(14kappa)<10^12/kappa.

The two proof regions cover every nonforward physical angle, including
the original gravity-dominated cone. Neither the matter amplitude nor
the Newton denominator has been expanded. Atkappa=10^800 the common
known rate bound is10^-788.

## Endpoint and unknown matching

For fixed original matter parameters, tau*A_m tends to zero and
tau*A_G tends toV/kappa. The leading real coefficient in notes/phase.md
therefore yields

I_known/sqrt(tau) ->3(5V+6)/(32kappa V).

The unknown finite Newton termdelta_kappa*T0/kappa^2 gives instead
a normalized interference limit-2delta_kappa/kappa. The local
alpha andbeta terms have zero normalized endpoint limit.
This does not bound or choose those physical coordinates.

The rate coefficient is in the fixed analytic hard reference.
Gravity-Born real radiation and its complete physical finite
conversion remain to be included. The Born cross section itself
is forward divergent. No loop square, all-order expansion or complex
Regge contour bound follows from the real one-loop estimate.
