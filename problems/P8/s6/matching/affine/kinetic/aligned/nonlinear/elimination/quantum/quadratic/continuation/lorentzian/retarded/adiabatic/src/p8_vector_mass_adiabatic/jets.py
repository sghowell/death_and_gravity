"""Exact commuting time jets at a fixed physical background metric."""
import sympy as sp
from p8_vector_dimensional import local

D, z = local.dimension, local.z
H = (local.H, local.Hd, local.Hdd, local.Hddd, local.Hdddd)+sp.symbols("H_fifth H_sixth", real=True)
alpha = (local.alpha,)+sp.symbols("alpha1:7", real=True)
beta = (local.beta,)+sp.symbols("beta1:7", real=True)
n = sp.symbols("source0:7", real=True)


def time(value):
    result = -2*H[0]*z*(1-z)*sp.diff(value, z)
    for row in (H, alpha, beta, n):
        result += sum(sp.diff(value, row[j])*row[j+1] for j in range(len(row)-1))
    return sp.expand(result)


def kind(value):
    if type(value) is not str or value not in ("T", "L"):
        raise ValueError("Require physical T or L sector")
    return value


def order(value):
    if type(value) is not int or value not in (0, 1, 2):
        raise ValueError("Require native adiabatic half-order 0,1 or 2")
    return value
