"""Order-zero exact overlap maps preserving the signed high-frequency subspaces."""
from functools import cache

import sympy as sp
from p8_proca_global_support import model, transfer

from . import canonical

clean=canonical.clean


@cache
def data():
    d=transfer.data()
    gamma=d["gamma"]["original_phase_to_scaled_velocity"]
    unitary=d["unitary"]["original_phase_to_scaled_velocity"]
    inverse_gamma=d["gamma"]["scaled_velocity_to_original_phase"]
    inverse_unitary=d["unitary"]["scaled_velocity_to_original_phase"]
    canonical_inverse=canonical.data()["packet_to_density"].subs(canonical.k,model.k)
    maps={"gamma_to_unitary":clean(unitary*inverse_gamma),
        "unitary_to_gamma":clean(gamma*inverse_unitary),
        "canonical_packet_to_gamma_velocity":clean(gamma*canonical_inverse)}
    out={}
    for name,T in maps.items():
        Ti=clean(T.inv())
        out[name]={"complete_map":T,"complete_inverse":Ti,
            "leading_map":clean(T.applyfunc(lambda v:sp.limit(v,model.k,sp.oo))),
            "leading_inverse":clean(Ti.applyfunc(lambda v:sp.limit(v,model.k,sp.oo))),
            "maximum_frequency_degree":max(transfer.rational_degree(v,model.k) for v in T),
            "inverse_maximum_frequency_degree":max(transfer.rational_degree(v,model.k) for v in Ti)}
    return out


@cache
def checks():
    d=data()
    o=model.old
    W=sp.diag(-o.theta/o.lam,1,-o.theta/o.lam,1)
    A=canonical.data()["leading_generator"][:2,2:]
    expected=sp.diag(*(o.a**sp.Rational(-3,2) for _ in range(4)))
    expected[2:,2:]=o.a**sp.Rational(-3,2)*A
    out={
        "global_gamma_unitary_leading_overlap_is_signed_frequency_preserving":clean(
            d["gamma_to_unitary"]["leading_map"]-W),
        "global_canonical_gamma_leading_velocity_map_is_positive_kinetic":clean(
            d["canonical_packet_to_gamma_velocity"]["leading_map"]-expected),
        "global_gamma_unitary_overlap_is_the_exact_inverse_pair":clean(
            d["gamma_to_unitary"]["complete_map"]*d["unitary_to_gamma"]["complete_map"]-sp.eye(4)),
    }
    for name,row in d.items():
        out[name+"_leading_overlap_is_invertible"]=clean(row["leading_map"]*row["leading_inverse"]-sp.eye(4))
    return out


@cache
def gates():
    return {name+"_and_inverse_have_order_zero_classical_symbols":bool(
        row["maximum_frequency_degree"]<=0 and row["inverse_maximum_frequency_degree"]<=0)
        for name,row in data().items()}
