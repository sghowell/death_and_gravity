"""Exact univariate derivatives encoded by symmetric labelled nilpotent jets."""
import sympy as sp
from p8_physical import jets as j
from sympy.polys.domains import QQ_I


def context(parameter,degree=4):
    if type(degree) is not int or not 1<=degree<=4:
        raise ValueError("Require native Taylor degree one through four")
    ctx=j.Context(((0,0,0),)*degree)
    ctx.domain=QQ_I.frac_field(parameter)
    ctx.zero_coefficient=ctx.domain.zero
    ctx.one_coefficient=ctx.domain.one
    return ctx,sum((ctx.leg(n) for n in range(degree)),ctx.jet())


def power(value,exponent):
    exponent=sp.sympify(exponent)
    if not isinstance(exponent,sp.Rational):
        raise TypeError("Require an exact rational exponent")
    constant=value.coefficient(0)
    if constant==0:
        if not isinstance(exponent,sp.Integer) or exponent<0:
            raise ValueError("Nonanalytic power at zero")
        return value**int(exponent)
    normalized=value/constant
    return constant**exponent*normalized.power(exponent)


def evaluate(expression,mapping,ctx):
    """No numerical sampling: recursively evaluate rational/algebraic expressions."""
    memo={}
    def visit(expr):
        if expr in memo:
            return memo[expr]
        if expr in mapping:
            value=mapping[expr]
            result=value if isinstance(value,j.Jet) else ctx.jet(value)
        elif not expr.free_symbols.intersection(mapping):
            result=ctx.jet(expr)
        elif expr.is_Add:
            result=sum((visit(arg) for arg in expr.args),ctx.jet())
        elif expr.is_Mul:
            result=ctx.jet(1)
            for arg in expr.args:
                result*=visit(arg)
        elif expr.is_Pow:
            result=power(visit(expr.base),expr.exp)
        else:
            raise TypeError("Unsupported Taylor expression: "+str(expr))
        memo[expr]=result
        return result
    return visit(sp.sympify(expression))


def derivatives(value):
    return tuple(sp.factor(value.coefficient((1<<degree)-1))
                 for degree in range(value.context.n+1))


def from_derivatives(ctx,epsilon,values):
    if len(values)>ctx.n+1:
        raise ValueError("Too many derivatives for this labelled Taylor context")
    return sum((ctx.jet(value)*epsilon**n/sp.factorial(n) for n,value in enumerate(values)),ctx.jet())


def integrate_zero(value,epsilon):
    """Primitive in epsilon, with exactly zero fixed basepoint."""
    values=derivatives(value)
    return from_derivatives(value.context,epsilon,(0,)+values[:-1])
