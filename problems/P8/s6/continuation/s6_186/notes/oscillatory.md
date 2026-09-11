# Full pair smearing and internal momentum integrals

Take any real smooth compact symmetric physical-frame
test tensor f_ab(t,x) in the open unit CD slab. Fourier
transform only space, with fhat(P)=integral e^(-iPx)f(x)dx
and integration measure d^3P/(2pi)^3.

For each pair let P=k+l, nu=sqrt(m^2+k^2/Amax^2) and
mu=sqrt(m^2+l^2/Amax^2). The pure-reference creation
pair has phase exp(i integral(W_k+W_l)), vector amplitude
norm at most C sqrt(nu mu), C=1e9, and inverse phase
g=1/(W_k+W_l) bounded by G=2/(nu+mu).
Both are analytic on radius r=1/20000 discs.

Three integrations by parts have no endpoint terms.
Ignoring unit-modulus powers of i, the adjoint operator
is L(b)=(g b)'. Its complete third iterate is

    L^3 b=g^3 b'''+6g^2 g' b''
      +(7g(g')^2+4g^2 g'')b'
      +((g')^3+4g g' g''+g^2 g''')b.

Apply this to amplitude times fhat. Banach-valued
Cauchy inequalities bound the amplitude and inverse
phase derivatives by j! C/r^j and j! G/r^j.
The coefficients of test derivatives of orders0..3
are bounded by(48,33,9,1) C G^3/r^(3-j).
Their sum is at most91 C G^3/r^3.

Cauchy-Schwarz on an interval of length1 consequently
gives, with Uj(P)^2=sum_(r<=j) integral
||partial_t^r fhat(t,P)||F^2 dt,

    |I_ref(k,l)|<=1e25 sqrt(nu mu)/(nu+mu)^3 U3(P).

The complete non-reference pair instead obeys

    |I_err(k,l)|<=1e15 sqrt(nu mu)
                   (nu^-6+mu^-6) U0(P).

No integration by parts is applied to its mixing
coefficients. This is not a cutoff in time or momentum.

For fixed P, the reference squared internal weight
satisfies nu mu/(nu+mu)^6<=1/(4nu^4).
The exact three-dimensional integrals are

    J4=integral nu^-4 d^3k/(2pi)^3=Amax^3/(8pi m),
    J10=integral nu^-10 d^3k/(2pi)^3
        =5Amax^3/(512pi m^7).

They follow from radial integrals pi/4 and5pi/256.
In particular J4/4<1/m and J10<1/m^7.

The Euclidean norm triangle inequality gives
mu<=nu+|P|/Amax<=nu L(P), with
L(P)=1+|P|/(Amax m). Symmetrizing the complete squared
remainder weight then gives an internal integral
at most4L(P)J10<8m^-7(1+|P|^2).
This uses the full k space and requires only a
spatial H1 test weight. It does not assume P=0.
