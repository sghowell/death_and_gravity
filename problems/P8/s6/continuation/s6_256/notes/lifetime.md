# Evaluated unforced classical lifetime with the full heavy mass

Use S6.255's entire normalized homogeneous action per a_hat^3,

    L=-3 D Hhat^2+3 B Hhat
      +N U F+9 U R_u^2/(16 R N)-I_u
      +Z(mc^2+mh^2)/2+N U[-n hbar^2/2+j hbar],
    B=-U R_u/(2N)-I,  I_N=3 U R_u R_N/(4 R N).

Here j is the COMPLETE normalized source, not just its bare term.
The code re-enters this literal action independently and removes the
primitive only by its exact lapse derivative. Undifferentiated I and
I_u cancel from the volume Euler expression with its a_hat^3 weight.

Let v=(Hhat,mc,mh), K0=diag(-6D,Z,Z), E=L_N and

    Theta=-Hhat D_N+B_N/2,
    q=E_v=(6Theta,Z_N mc,Z_N mh),
    Cnn=E_N/2,
    J=Cnn-q^T K0^-1 q/2,
    f=(3L-3Hhat L_Hhat-partial_u L_Hhat,
       -3Hhat Z mc-Z_u mc,
       N U(-n hbar+j)-3Hhat Z mh-Z_u mh),
    rest=E_u+E_hbar mh.

The complete Euler/constraint system is

    K0 v_dot+q N_dot=f,
    q^T v_dot+2 Cnn N_dot=-rest.

Its exact solution is

    N_dot=[-rest-q^T K0^-1 f]/(2J),
    v_dot=K0^-1(f-q N_dot).

Every component is retained. The action Hessian and cross vector are
checked directly; an independent symbolic four-by-four Schur identity
checks the entire solve. The initial constraint is zero, and the
differentiated constraint remains zero. The unchanged S6.255 reconstruction
then gives the lapse, volume, both matter, vector, affine and clock
Euler equations. This is the local classical action, not an assertion
that its coefficients are the self-consistent interacting quantum mean.

Bootstrap on the connected solution component containing the initial
datum, with |u|<=.001, |N-1|<=.001, .999<=a_hat<=1.001,
.99<=R<=1.01, 1<=J<=2, |Hhat|<=.001, .09<=mc<=.11,

    |mh|<=10^-1000, |hbar|<=10^-1098,
    10^196<=n<=10^198,
    all full R four-jets <=200, F four-jets <=10^4,
    all full j four-jets <=10^-2500.

Exact outward enclosures of the ENTIRE equations give
|N_dot|,|Hhat_dot|,|mc_dot|,|mh_dot|<10^20.
The full first time derivatives of
D,Z,J,Cnn,Theta,r,R_N,Y,K,Yv,L2,C,c_M1,w_M1
are below 10^40. Derivatives include every state rate, explicit source
jet, lapse dependence, and a_hat_dot=a_hat Hhat. Profiles remain fixed
functions of time but their derivatives are not erased.

A raw mass-squared Lipschitz bound would be a useless estimate here.
Instead retain the positive heavy energy

    e^2=mh^2+N^2 n hbar^2,
    mh_dot=-(3Hhat+Z_dot/Z)mh-N^2 n hbar+N^2 j.

Its exact derivative is

    (e^2)'=-2(3Hhat+Z_dot/Z)mh^2
            +2N^2 j mh+2N N_dot n hbar^2.

The mass cross terms cancel, not the mass or source. The actual
rational N_dot,Z_dot enclosures give root-energy damping below 10^20.
Since N^2<2, the upper derivative of e satisfies

    e' <= 10^20 e+2*10^-2500,  e(0)=0.

The inequality at e=0 follows by regularizing the energy or the usual
upper-Dini-derivative argument. For T=10^-60, exp(10^20 T)<2, hence
e<=4 T*10^-2500=4*10^-2560. This strictly improves both heavy bootstrap
bounds: the actual n>10^197 and N>.999 imply N sqrt(n)>10^98.

Light rate and lapse displacements are below 10^-40; all persistent
coefficient displacements are below 10^-20. The 14 computed initial
margins exceed twice that latter bound. The scale factor remains
strictly inside its box by a_hat_dot=a_hat Hhat. The clock interval
lies strictly inside the full coefficient strip. Theta starts at zero,
so |Theta|<10^-20; |Hhat-Hclock|<10^-40+4T<10^-20, where
Hclock=4u/(1+u^2).

The improved energy, together with full coefficient enclosures, gives
|c_H|,|w_H|<10^-1000. The ENTIRE mixed potential source
partial_N[N U(-n hbar+j)] is enclosed below 10^-880 on the bootstrap
box. Also mu=N U n stays in [n/2,2n] and |mu_dot/mu|<10^40.
All these are checked against the exact matrix input box.

If the solution stopped before T, these strict improvements would put
it in a compact subset of the regular J,Cnn,positive-coefficient domain.
The full finite-dimensional vector field and its needed derivatives
remain bounded there. Local ODE continuation extends the solution,
contradicting the alleged endpoint. Thus the full-current unforced
classical comparison reaches T. This is an evaluated interval at the
fixed giant mass, not a uniform-in-mass theorem or numerical integration
of that physical oscillator, and not a macroscopic physical bounce.
