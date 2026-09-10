# Uniform isolated coordinate bound

For M>=1 split J1 at y=M. On 0<=y<=M, 1/(y+M)<=1/M.
On y>=M, its full integrand is <=1/y^2. Integrating gives

    J1 <= [log(M+1)+1/(M+1)]/M <= log(4M)/M.

For the second inequality use M+1<=2M, 1/(M+1)<=1/2 and
log 2=integral_1^2 dt/t>=1/2. At the actual parameters
4M<10^200 and log 10<3 (exp 3>1+3+9/2+27/6>10),
so log(4M)<600. Together with g/M<L/3 and Q>144 this yields

    0<k<25L/3.

The already-established fixed-contact bound is
|sigma|<=50L^2/3. Consequently, for 0<=h<=1 and
0<L<3/50, 1-hk>1/2. The inverse is nonsingular, and

    |h^2 sigma partial_L sigma| <= h^2 (1250/9)L^3,
    |Lstar-L+h sigma-h^2 sigma partial_L sigma|
       <= h^3 |sigma| k^2/(1-hk).

The exact signed remainder is -h^3 sigma k^2/(1-hk).
At the actual parameters the bounds for the first shift,
second shift and remaining coordinate tail are approximately
8.73114913702011e-409, 1.66533453693773e-612 and
3.17637355220363e-816, respectively. Since sigma<0,
Lstar-L=-h sigma/(1-hk)>=0: this isolated map cannot worsen
L-3g/M. The full combined map is NOT covered by that assertion.

Before matching, the assigned insertion has the S6.125 bound
(10025/18)L^3<7e-612. Its paired nonlocal order-two contribution
is exactly zero, not a triangle-inequality estimate. The inverse
coordinate tail is not a bound on physical higher-loop graphs,
nor may this identity alone subtract an item from a larger
budget whose other scheme-conversion terms are still missing.
