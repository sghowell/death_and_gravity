"""Independent pure-vector quartic anchor, reality and permutation tests."""
from functools import cache

import sympy as sp
from p8_auxiliary_neighborhood import model
from p8_coupled_momentum import constraints
from p8_physical import jets as j

from . import physical


@cache
def pure_electric():
    directions=((1,0,0),(0,1,0),(1,1,0),(-2,-2,0))
    ctx=constraints.context(directions)
    f={name:ctx.jet() for name in ("zeta","scalar_p","chi","chi_p")}
    t,tp=j.zeros(ctx),j.zeros(ctx)
    W=[ctx.jet() for _ in range(3)]
    Pi=[ctx.jet(),ctx.jet(),sum((ctx.leg(n) for n in range(4)),ctx.jet())]
    data=physical.construct(ctx,f,t,tp,W,Pi,point=0)
    Je=sp.Rational(1199,800)+4*model.MARGIN
    expected=sp.Rational(3,8)/(Je*sp.Rational(1,10**6)**2)
    return {"kernel":data["full_labelled_kernel"],"independent_kernel":expected,
            "omitting_lapse_response_quartic_defect":-expected,
            "checks":data["reduction"]["checks"],
            "residual":sp.factor(data["full_labelled_kernel"]-expected)}


def transformed_fixture(n,permutation=None,reverse=False):
    if type(n) is not int or n not in (3,4):
        raise ValueError("Require cubic or quartic anchor")
    if permutation is None:
        permutation=tuple(range(n))
    if type(permutation) is not tuple or any(type(x) is not int for x in permutation) or sorted(permutation)!=list(range(n)):
        raise ValueError("Require a native permutation of all external labels")
    if type(reverse) is not bool:
        raise TypeError("Require a native reality flag")
    return _transformed_fixture(n,permutation,reverse)


@cache
def _transformed_fixture(n,permutation,reverse):
    original=physical.fixture(n)
    old=original["context"]
    sign=-1 if reverse else 1
    directions=tuple(tuple(sign*x for x in old.momenta[index]) for index in permutation)
    ctx=constraints.context(directions)
    f,t,tp,W,Pi=physical.mixed_fields(old)
    def change(value):
        return sum((ctx.leg(new,value.coefficient(1<<before))
                     for new,before in enumerate(permutation)),ctx.jet())
    f={name:change(value) for name,value in f.items()}
    t,tp=([[change(value) for value in row] for row in matrix] for matrix in (t,tp))
    W,Pi=([change(value) for value in vector] for vector in (W,Pi))
    actual=physical.construct(ctx,f,t,tp,W,Pi,point=0)["full_labelled_kernel"]
    expected=sp.conjugate(original["full_labelled_kernel"]) if reverse else original["full_labelled_kernel"]
    return {"actual":actual,"expected":expected,"residual":sp.factor(actual-expected)}


@cache
def checks():
    out={"independent_pure_transverse_electric_quartic":pure_electric()["residual"]}
    for n in (3,4):
        out[f"mixed_{n}_leg_full_momentum_reversal_reality"]=transformed_fixture(n,reverse=True)["residual"]
        out[f"mixed_{n}_leg_all_label_reversal_permutation"]=transformed_fixture(n,permutation=tuple(reversed(range(n))))["residual"]
    return out
