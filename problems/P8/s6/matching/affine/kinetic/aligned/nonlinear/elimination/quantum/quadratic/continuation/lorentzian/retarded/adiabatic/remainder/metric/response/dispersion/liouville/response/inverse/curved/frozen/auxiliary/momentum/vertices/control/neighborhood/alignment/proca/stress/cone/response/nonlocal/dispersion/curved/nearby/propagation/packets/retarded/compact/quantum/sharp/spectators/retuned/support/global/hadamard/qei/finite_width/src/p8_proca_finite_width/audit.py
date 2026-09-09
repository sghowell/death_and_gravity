"""Independent exact identities and whole-interval analytic inputs."""
from functools import cache

from . import proofs


@cache
def residuals():
    return proofs.checks()

@cache
def gates():
    return {name:bool(value) for name,value in proofs.gates().items()}

@cache
def controls():
    return proofs.controls()
