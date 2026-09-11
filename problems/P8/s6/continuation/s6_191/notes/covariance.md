# Full real covariance derivatives

The actual and reference graphs are complex symmetric
and lie in the common radius1/10 ball. Treat their
entries as a real Banach space, so differentiating an
adjoint is included, not suppressed.

Set C(r)=(I-rdag r)^-1. At rho=1/10 its first three
Frechet derivative bounds follow from the ordered
inverse rule. Each derivative of rdag r of order1
costs2rho, of order2 costs2, and of higher order is
zero. With ||C||<2, the bounds are

    ||DC|| <=4(2rho)<1,
    ||D2C|| <=2*8(2rho)^2+4*2<16,
    ||D3C|| <=6*16(2rho)^3+6*8*2(2rho)<32.

The order and counts include both product directions.
For F(r)=[I+r; -i(I-r)]/sqrt(2), ||F||<2 and ||DF||=1,
while its higher derivatives vanish. The complete
six-quadrature map Sigma=Re(F C Fdag) consequently has

    ||D Sigma||<16, ||D2 Sigma||<128, ||D3 Sigma||<512.

For example the second derivative has one C'' term,
four F'/C' terms and two F'/F' terms. The third has
one C''' term, six F'/C'' terms and six F'/C'/F' terms.
The code checks these factors and the full inverse
and covariance jets on noncommuting complex matrices.

Write Delta=Sigma(r)-Sigma(rhat). The mean-value
identity on the convex graph ball and the ordinary
parameter chain rule give

    ||Delta||<=16E0,
    ||Delta1||<=16E1+128E0||rhat1||,
    ||Delta2||<=16E2+128E0||rhat2||
        +128E1(||r1||+||rhat1||)+512E0||rhat1||^2.

Use ||r1||<=||rhat1||+E1; the E1^2 term is retained.
The resulting exact constants, converted with the
common K lower bound only where necessary, imply

    ||Delta||<2e31 nu_minus^-10,
    ||Delta1||<4e33 nu_minus^-9,
    ||Delta2||<2e36 nu_minus^-8.

These are derivatives of the full balanced covariance
difference, including cross-correlations. They are
not derivatives of an independently reset Gaussian
state or a response of the gravitational metric itself.
