"""Physical H2,H3,H4 from the complete spatial and auxiliary reduction."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_momentum import constraints, quadratic
from p8_physical import geometry
from p8_physical import jets as j
from p8_s5.lapse_series import stationary_series

from . import coefficients


def background(point=None):
    if point is not None:
        point=constraints.exact(point)
        if abs(point)>sp.Rational(1,2):
            raise ValueError("Require rational time in I")
    return {key:value if point is None else sp.factor(value.subs(model.u,point))
            for key,value in model.coefficients()["background"].items()}


def invariants(ctx,geo,momentum,f,t,tp,W,Pi,bg):
    volume=geo["volume"]
    mixed=j.mscale(j.mmul(momentum,geo["metric"]),volume.power(-1))
    trace=j.trace(mixed)
    ell,H=bg["ell"],bg["H"]
    zeta=sp.Rational(1,10**6)
    F=[[W[b].derivative(a)-W[a].derivative(b) for b in range(3)] for a in range(3)]
    values=(
        sp.Rational(2,3)*trace+2*H,
        (ell+f["chi_p"])/volume-ell,
        sum(Pi[a].derivative(a) for a in range(3))/volume,
        j.trace(j.mmul(mixed,mixed))-trace*trace/3,
        sum(geo["metric"][a][b]*Pi[a]*Pi[b] for a in range(3) for b in range(3))/zeta/volume**2,
        zeta*sum(geo["inverse"][a][c]*geo["inverse"][b][d]*F[a][b]*F[c][d]
                  for a in range(3) for b in range(3) for c in range(3) for d in range(3)),
        sum(geo["inverse"][a][b]*W[a]*W[b] for a in range(3) for b in range(3)),
        sum(geo["inverse"][a][b]*f["chi"].derivative(a)*f["chi"].derivative(b)
              for a in range(3) for b in range(3)),
        geo["curvature"])
    if any(not value.homogeneous(0).is_zero() for value in values):
        raise ValueError("A physical invariant has nonzero background")
    if any(not values[n].homogeneous(1).is_zero() for n in range(3,8)):
        raise ValueError("A squared physical invariant has a linear term")
    return values


def invariant_reduction(ctx,values,point=None):
    rows=coefficients.jets()
    zero=(0,)*9
    evaluate=lambda x:x if point is None else sp.factor(x.subs(model.u,point))
    A=tuple(evaluate(value) for value in rows[zero])
    if A[1]!=0 or A[2]==0:
        raise ValueError("Require the actual stationary invertible lapse background")
    L=[ctx.jet() for _ in range(4)]
    Q=[ctx.jet() for _ in range(3)]
    for powers,row in rows.items():
        if sum(powers)==0:
            continue
        if sum(powers) not in (1,2):
            raise ValueError("Unexpected invariant polynomial order")
        term=ctx.jet(1)
        for value,power in zip(values,powers):
            term*=value**power
        destination=L if sum(powers)==1 else Q
        for degree in range(len(destination)):
            destination[degree]+=evaluate(row[degree])*term
    if ctx.n==2:
        # A quadratic Hamiltonian needs only the linear physical lapse.
        # Avoid forming higher rational lapse coefficients that cannot
        # enter its two-label algebra.
        n1=-L[1].homogeneous(1)/A[2]
        result={"lapse":(n1,),"hamiltonian":(A[0],L[0],Q[0]-L[1]*L[1]/(2*A[2]))}
    else:
        result=stationary_series(A,L,Q)
    h=sum(result["hamiltonian"],ctx.jet())
    # This checks the actual physical stationary force, not only a formal
    # invariant-weight fixture. Higher physical invariant terms are retained.
    n=sum(result["lapse"],ctx.jet())
    force=A[2]*n+A[3]*n*n/2+A[4]*n**3/6+L[1]+L[2]*n+L[3]*n*n/2+Q[1]+Q[2]*n
    direct=A[0]+A[2]*n*n/2+A[3]*n**3/6+A[4]*n**4/24
    direct+=L[0]+L[1]*n+L[2]*n*n/2+L[3]*n**3/6
    direct+=Q[0]+Q[1]*n+Q[2]*n*n/2
    checks={f"actual_physical_lapse_force_degree_{degree}":force.homogeneous(degree).is_zero()
            for degree in range(ctx.n)}
    checks["direct_stationary_Hamiltonian_through_four"]=(direct-h).is_zero()
    if not all(checks.values()):
        raise ValueError("Actual physical lapse reconstruction failed")
    return {"density":h,"series":result,"checks":checks,"A":A,"L":tuple(L),"Q":tuple(Q),
            "stationary_force":force,"direct_Hamiltonian_difference":direct-h}


def construct(ctx,f,t,tp,W,Pi,point=None):
    fields=tuple(f.values())+tuple(x for matrix in (t,tp) for row in matrix for x in row)+tuple(W)+tuple(Pi)
    if any(value.context is not ctx or not value.homogeneous(0).is_zero()
           or not (value-value.homogeneous(1)).is_zero() for value in fields):
        raise ValueError("Require zero-background linear fields in one labelled context")
    for matrix in (t,tp):
        if not j.trace(matrix).is_zero() or any(not (matrix[a][b]-matrix[b][a]).is_zero()
                                                for a in range(3) for b in range(3)):
            raise ValueError("Require symmetric tracefree tensor fields and momenta")
        if any(not sum(matrix[a][b].derivative(b) for b in range(3)).is_zero() for a in range(3)):
            raise ValueError("Require transverse tensor fields and momenta")
    bg=background(point)
    geo=geometry.derive(ctx,f["zeta"],t)
    solved=constraints.derive(ctx,geo,f["zeta"],t,f["scalar_p"],tp,
        f["chi"],f["chi_p"],W,Pi,bg["H"],bg["ell"])
    values=invariants(ctx,geo,solved["momentum"],f,t,tp,W,Pi,bg)
    reduced=invariant_reduction(ctx,values,point)
    H=reduced["density"]*geo["volume"]-2*bg["H"]*j.contract(solved["momentum"],geo["metric"])
    Adot=sp.diff(model.coefficients()["background"]["H"],model.u)+3*model.coefficients()["background"]["H"]**2
    if point is not None:
        Adot=sp.factor(Adot.subs(model.u,point))
    H-=Adot*(6*f["zeta"]+3*f["zeta"]**2-j.contract(t,t)/2)+bg["ell"]*f["chi_p"]
    return {"context":ctx,"background":bg,"geometry":geo,"momentum":solved,
            "invariants":values,"reduction":reduced,"Hamiltonian":H,
            "full_labelled_kernel":sp.factor(H.coefficient(ctx.full))}


def mixed_fields(ctx):
    f={name:ctx.jet() for name in ("zeta","scalar_p","chi","chi_p")}
    t,tp=j.zeros(ctx),j.zeros(ctx)
    W,Pi=([ctx.jet() for _ in range(3)] for _ in range(2))
    for leg,wave in enumerate(ctx.momenta):
        k=sp.Matrix(wave)
        A=next(k.cross(sp.eye(3)[:,a]) for a in range(3) if k.cross(sp.eye(3)[:,a])!=sp.zeros(3,1))
        B=k.cross(A)
        E=A*A.T/(A.T*A)[0]-B*B.T/(B.T*B)[0]
        for index,name in enumerate(f):
            f[name]+=ctx.leg(leg,sp.Rational(leg+index+1,17+2*index))
        for a in range(3):
            W[a]+=ctx.leg(leg,sp.Rational(leg+a+2,29))
            Pi[a]+=ctx.leg(leg,sp.Rational((leg+1)*(a+1)+1,37))
            for b in range(3):
                t[a][b]+=ctx.leg(leg,E[a,b]*sp.Rational(leg+1,31))
                tp[a][b]+=ctx.leg(leg,E[a,b]*sp.Rational(leg+1,82))
    return f,t,tp,W,Pi


def fixture(n,point=0,*,symbolic_scale=False):
    if type(n) is not int or n not in (2,3,4):
        raise ValueError("Require native physical degree two, three or four")
    if type(symbolic_scale) is not bool:
        raise TypeError("Require a native scale flag")
    if point is not None:
        point=constraints.exact(point)
    background(point)
    return _fixture(n,point,symbolic_scale)


@cache
def _fixture(n,point,symbolic_scale):
    directions={2:((1,2,1),(-1,-2,-1)),
                3:((1,0,0),(0,1,0),(-1,-1,0)),
                4:((1,0,0),(0,1,0),(0,0,1),(-1,-1,-1))}[n]
    ctx=constraints.context(directions,symbolic_time=point is None,symbolic_scale=symbolic_scale)
    return construct(ctx,*mixed_fields(ctx),point)


def quadratic_pair(left,right,point=0):
    if type(left) is not str or type(right) is not str or left not in quadratic.CHANNELS or right not in quadratic.CHANNELS:
        raise ValueError("Require declared quadratic phase channels")
    point=constraints.exact(point)
    background(point)
    return _quadratic_pair(left,right,point)


@cache
def _quadratic_pair(left,right,point):
    ctx=constraints.context(((1,0,0),(-1,0,0)),symbolic_scale=True)
    data=construct(ctx,*quadratic.fields(ctx,left,right),point=point)
    return sp.factor(data["full_labelled_kernel"]-quadratic.pair(left,right)["actual"].subs(model.u,point))


@cache
def quadratic_algebraic_bridge():
    # Equality at invariant level proves the quadratic identity for all
    # modes and directions, before the common spatial substitution.
    bg=model.coefficients()["background"]
    rows=coefficients.jets()
    Fclock=sp.factor(model.adm.coefficients()["F"].subs(model.adm.X,1))
    expected=(-sp.Rational(3,4)*(model.dp-2*bg["H"])**2-Fclock+model.j**2/2
              +2*model.shear+model.electric/2+model.magnetic/4+model.vector/2
              +(bg["ell"]+model.dc)**2/2+model.chi_gradient/2-model.curvature/2)
    point=dict(sp.Poly(sp.expand(expected),*model.VARIABLES).terms())
    out={}
    for powers in set(rows)|set(point):
        out["full_quadratic_on_clock_invariant_"+"_".join(map(str,powers))]=sp.factor(
            rows.get(powers,(sp.Integer(0),)*5)[0]-point.get(powers,sp.Integer(0)))
    first={model.dp:3*bg["theta"],model.dc:-bg["w"],model.j:0,model.curvature:-bg["lam"]/2}
    for variable,value in first.items():
        powers=tuple(int(v==variable) for v in model.VARIABLES)
        out["full_quadratic_linear_lapse_"+str(variable)]=sp.factor(rows[powers][1]-value)
    zero=(0,)*9
    out["full_quadratic_stationary_background"]=rows[zero][1]
    out["full_quadratic_lapse_Hessian"]=sp.factor(rows[zero][2]+2*quadratic.compact_scalar()["J_e"])
    return out
