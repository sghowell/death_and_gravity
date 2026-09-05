"""Three-dimensional metric geometry for the exact linear spatial gauge."""

from . import jets as j


def derive(context, zeta, tensor):
    unit = j.identity(context)
    perturbation = j.madd(j.mscale(unit, 2*zeta), tensor)
    metric = j.madd(unit, perturbation)
    inverse = unit
    term = unit
    for degree in range(1, context.n+1):
        term = j.mmul(term, perturbation)
        inverse = j.madd(inverse, j.mscale(term, (-1)**degree))
    volume = j.determinant(metric).power(j.sp.Rational(1, 2))
    # Gamma^i_jk, constructed from the coordinate metric rather than an R shortcut.
    lower = [[[perturbation[a][k].derivative(i) for k in range(3)] for i in range(3)] for a in range(3)]
    christoffel = [[[sum(inverse[i][l]*(lower[l][k][a]+lower[l][a][k]-lower[a][l][k])/2
                        for l in range(3)) for k in range(3)] for a in range(3)] for i in range(3)]
    # lower[l][k][a]=partial_k g_la; lower[a][l][k]=partial_l g_ak.
    ricci = [[sum(christoffel[k][a][b].derivative(k)-christoffel[k][a][k].derivative(b)
                  +sum(christoffel[k][a][b]*christoffel[l][k][l]
                       -christoffel[l][a][k]*christoffel[k][b][l] for l in range(3))
                  for k in range(3)) for b in range(3)] for a in range(3)]
    curvature = j.contract(inverse, ricci)
    return {"metric": metric, "inverse": inverse, "volume": volume,
            "christoffel": christoffel, "curvature": curvature}
