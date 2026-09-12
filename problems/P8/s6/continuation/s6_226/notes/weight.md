# Joint time-space frequency estimate and the original causal multiplier

Let lambda=sigma+i omega, sigma>0, q>=0 andp=lambda^2+q. Direct algebra gives

|p|^2-sigma^2(sigma^2+omega^2+q)
=(q-omega^2)^2+sigma^2(q+omega^2)>=0.

Consequently |p|>=sigma^2 even near the wave coneq=omega^2. If omega is nonzero, Im p=2sigma omega is nonzero; if omega=0, p=sigma^2+q>0. The whole frequency line stays on the original uncut first sheet.

For M>=0 choose sigma=2m exp(32+8M). Then

log(sigma^2/(4m^2))=64+16M,
Re A_i(lambda^2+q)>=(32+13M)/5=d_M>0.

The actual diagonal inverse multiplier is

diag(1/Ftrace(lambda^2+q),(3/8)/F2(lambda^2+q)),

so its Euclidean norm is at most1/d_M=5/(32+13M), uniformly in BOTH omega and spatial momentum. The factor3/8 is retained; no Frobenius or channel normalization is changed.

## From the original causal distribution to weighted L2

Use the already constructed original causal inverse kernels, not a newly chosen Green function. For every fixedq, K_i,q=D J_i,q with the S224 uniformly bounded primitive. Thus the causal distribution has a Fourier-Laplace transform for sigma>0. For smooth compactly time-supported and compact-momentum inputs, multiplying by exp(-sigma t), Fourier transforming time and space, and using the known reciprocal transform gives the displayed multiplier.

Plancherel proves the uniform weighted L2_tH^r estimate for every realr. Density extends it to the whole weighted L2 space. The subspace of functions vanishing before a given time is closed in that norm; approximation and the original causal support therefore preserve causality under the extension. The resulting bounded operator is the same original distributional convolution.

No uniform unweighted half-line L1 bound is required, and no shear pole is removed. The exponential factor is a norm used to control an already fixed causal operator, not a modification of its physical equation or a damping term inserted into the action.
