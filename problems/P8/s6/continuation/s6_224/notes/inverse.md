# Complete causal matrix inverse and exact convolution constants

Let C=5/2, Jmax=45/2. The coordinate inverse kernel Rq has Rq(0)=0, ||Rq'||<C and ||Rq(t)||<Ct. The diagonal primitive Jq is continuous with Jq(0)=0 and norm at mostJmax. Consequently

Eq=Rq'*Jq*Rq^T

is an ordinary matrix-valued causal convolution. No frequency-absolute scalar kernel integral or uniform scalarL1 estimate is used. The transpose acts only on matrix entries, not time orientation. All factors here are constant-coefficient flat reference convolution operators.

For Re(lambda)>0, L[Rq]=B(lambda)^-1 and L[Rq']=lambda B^-1. The primitive formula therefore gives

L[Eq]=B^-1 diag(1/Ftrace(lambda^2+q),(3/8)/F2(lambda^2+q)) (B^T)^-1.

The lambda factors cancel exactly. Both products with A_q=B^T diag(Ftrace,(8/3)F2)B equal the identity matrix. This is checked in both orders. The determinant is(8/3)Ftrace F2 lambda^4(lambda^2+q)^2, nonzero throughout Re(lambda)>0.

The forward scalar distributions were constructed explicitly, and B is a differential matrix. All relevant causal convolutions exist; their Laplace transforms are of at most polynomial growth on each fixed right half-plane. Transform uniqueness proves both A_q*Eq=I delta and Eq*A_q=I delta. The three-source gauge block has not been inverted.

The kernel norm follows from the actual ordered triangle integral:

||Eq(t)||<=C^2 Jmax integral_(u=0)^t integral_(v=0)^(t-u)(t-u-v)dvdu
         =375t^3/16.

Integrating in time gives the operator bound

integral_0^T ||Eq(t)||dt<=375T^4/64<6T^4.

These constants are uniform in every spatial momentum and the positive fixed mass. They are normalized reference bounds, not physical-force smallness.

There is also a uniform first-time-derivative bound. Differentiate the last Rq factor in the convolution; its vanishing initial value removes the endpoint term. No derivative of Jq or second derivative of Rq is required:

Eq'=Rq'*Jq*(Rq')^T,
||Eq'(t)||<=C^2 Jmax t^2/2=1125t^2/16,
integral_0^T ||Eq'(t)||dt<=375T^3/16.

The expression is continuous for eachq, and the same dominating constants give the functional-space extension. Higher time/spatial derivatives are not silently inferred from this first-derivative argument.

For a finite window, the causal graph consists of continuous quotient amplitudes for which the complete forward distribution equals an ordinary continuous dual source, including the initial boundary. Extending histories beyondT does not affect the equation beforeT. For every such source f, y=Eq*f lies in that graph and solves the equation; conversely Eq*(A_q*y)=y proves uniqueness. Distributional initial atoms are not discarded. The result concerns this graph, not density of an unrelated domain or graph invariance of the curved S222 force system.
