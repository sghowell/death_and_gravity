# Primary special-function and internal-source audit

The science here is a classical ODE connection and a finite-window transfer
bound. No literature cutoff, vacuum, particle-production interpretation, or
health verdict is imported.

1. [NIST DLMF 15.10.1](https://dlmf.nist.gov/15.10.E1) defines the ordinary
   Gauss equation. Substitution of the two parity ansatzes gives the exact
   numerator sums/products in `connection.identities()`. The argument
   \(-8x^2\) stays outside the principal Gauss branch cut. The right Jost
   ansatz instead uses parameters \((-1/2,3/2;1-i\rho)\) and argument
   \((1-\tanh t)/2\) strictly between zero and one.

2. [DLMF 15.10.21](https://dlmf.nist.gov/15.10.E21), together with the
   definitions 15.10.13--14, is the **ordinary** Gauss connection used for
   A and B. Its coefficients are
   \(\Gamma(c)\Gamma(c-a-b)/[\Gamma(c-a)\Gamma(c-b)]\) and
   \(\Gamma(c)\Gamma(a+b-c)/[\Gamma(a)\Gamma(b)]\). Here
   \(a+b=1\), \(c=1-i\rho\); the relevant 0-to-1 exponent difference
   \(c-a-b=-i\rho\) is noninteger and all displayed Gamma coefficients
   are finite. Although \(a-b=-2\) is integer at the unused infinity
   singularity, this connection extends to it by analytic continuation in
   the parameters. As \(t\to-\infty\), \(1-w\sim e^{2t}\), yielding
   exactly the time-frequency convention printed in the proof.

3. [DLMF 15.8.2](https://dlmf.nist.gov/15.8.E2) supplies an independent
   infinity connection for the parity basis. **Its bold F is Olver's
   normalized function, not ordinary \({}_2F_1\).** Converting
   \(\mathbf F(a,b;c;z)={}_2F_1(a,b;c;z)/\Gamma(c)\), including the changed
   lower parameters of the functions at \(1/z\), produces E and O in the
   proof. The relevant \(-z=8x^2\) is positive, so the stated principal
   branch condition is satisfied. Numerical parity/Jost agreement is an
   additional regression, not a replacement for this normalization audit.

4. [DLMF 5.5.3](https://dlmf.nist.gov/5.5.E3) gives Gamma reflection;
   [5.5.1](https://dlmf.nist.gov/5.5.E1) gives recurrence. In particular
   \(\Gamma(-1/2)\Gamma(3/2)=-\pi\), so B has the positive imaginary sign
   \(i\operatorname{csch}(\pi\rho)\) in the declared convention.
   [5.4.3--4](https://dlmf.nist.gov/5.4.E3) give the imaginary-axis and
   half-integer modulus identities. Applying recurrence to both denominator
   factors gives \(|A|^2=\coth^2(\pi\rho)\); the recurrence factors cancel.

The frozen S6.21 `p8_variable_response.connection` and its proof were read
for the generic phase/orientation dictionary. No import of its coupled-mode
\(\mu=\sqrt{39}/2\), its small-mixing bound, source functional, or transfer
verdict is made. The own-f coefficient is 16 and \(\rho=\sqrt7/2\).

The original physical-time and tensor-action normalization comes from the
frozen S6.33 stationary branch, not the special-function source. Within this
child, `potential.py` supplies the full actual canonical remainder bound;
`connection.py` is explicitly conditional on that premise until the wrapper
replays it. Neither the frozen branch nor this note identifies a fictitious
free exterior with a physical extension of the parent action.
