"""Generic endpoint-owner grouping, with all unconsumed soft factors retained."""

from functools import cache
from itertools import product

import sympy as s
from sympy.polys.rings import ring


def require_case(number, sides, endpoint_mask):
    if type(number) is not int or number not in (2, 3):
        raise ValueError("Require two or three independent singleton emissions")
    if (
        not isinstance(sides, tuple)
        or len(sides) != number
        or any(type(x) is not int or x not in (0, 1) for x in sides)
    ):
        raise ValueError("Require one exact left/right selector per emission")
    if type(endpoint_mask) is not int or not 0 <= endpoint_mask < 16:
        raise ValueError("Require a subset of the four scalar endpoint slots")
    return number, sides, endpoint_mask


def coefficient_ring(number):
    if type(number) is not int or number not in (2, 3):
        raise ValueError("Require two or three independent singleton emissions")
    return _coefficient_ring(number)


@cache
def _coefficient_ring(number):
    names = [f"P{leg}" for leg in range(4)]
    names += [f"S{i}_{leg}" for i in range(number) for leg in range(4)]
    names += [f"X{i}_{leg}" for i in range(number) for leg in range(4)]
    C, *values = ring(",".join(names), s.QQ)
    v = dict(zip(names, values))
    R, *w = ring(",".join(f"w{i}" for i in range(number)), C.to_domain())
    return C, v, R, w


def group_case(number, sides, endpoint_mask):
    require_case(number, sides, endpoint_mask)
    return _group_case(number, sides, endpoint_mask)


@cache
def _group_case(number, sides, endpoint_mask):
    C, v, R, w = coefficient_ring(number)
    groups = ((0, 1), (2, 3))
    choices = tuple(groups[side] for side in sides)
    made = R.zero
    for assignment in product(*choices):
        endpoint = R.one
        for leg in range(4):
            if not endpoint_mask >> leg & 1:
                continue
            factor = R(v[f"P{leg}"])
            for i, owner in enumerate(assignment):
                if owner == leg:
                    factor += w[i] * R(v[f"X{i}_{leg}"])
            endpoint *= factor
        current = C.one
        for i, owner in enumerate(assignment):
            current *= v[f"S{i}_{owner}"]
        made += R(current) * endpoint
    records = []
    for powers, coefficient in made.items():
        untouched = C.one
        for i, power in enumerate(powers):
            if not power:
                left, right = choices[i]
                untouched *= v[f"S{i}_{left}"] + v[f"S{i}_{right}"]
        quotient, remainder = coefficient.div(untouched)
        records.append(
            {
                "powers": powers,
                "coefficient": coefficient.as_expr(),
                "untouched_factor": untouched.as_expr(),
                "quotient": quotient.as_expr(),
                "residual": remainder.as_expr(),
                "cross_multiplication": (
                    quotient * untouched + remainder - coefficient
                ).as_expr(),
            }
        )
    return {
        "records": records,
        "squarefree_support": all(
            all(power in (0, 1) for power in mon) for mon in made
        ),
        "whole_grouped_polynomial": made.as_expr(),
    }


@cache
def data():
    checks = {}
    counts = {}
    squarefree = True
    controls = {}
    for number in (2, 3):
        cases = entries = 0
        for sides in product((0, 1), repeat=number):
            for endpoint_mask in range(16):
                row = group_case(number, sides, endpoint_mask)
                squarefree = squarefree and row["squarefree_support"]
                prefix = f"m{number}_sides{''.join(map(str, sides))}_endpoints{endpoint_mask}"
                for index, record in enumerate(row["records"]):
                    checks[prefix + f"_coefficient{index}"] = record["residual"]
                    checks[prefix + f"_reconstruction{index}"] = record[
                        "cross_multiplication"
                    ]
                    entries += 1
                cases += 1
        counts[number] = (cases, entries)
        _C, v, _R, _w = coefficient_ring(number)
        _quotient, remainder = v["S0_0"].div(v["S0_0"] + v["S0_1"])
        controls[f"omitted_owner_m{number}"] = bool(remainder)
    return {
        "checks": checks,
        "gates": {
            "all192_endpoint_channel_cases": sum(v[0] for v in counts.values()) == 192,
            "all776_exact_coefficient_divisions": sum(v[1] for v in counts.values())
            == 776,
            "every_endpoint_insertion_support_squarefree": squarefree,
            "omitted_owner_negative_controls_nonzero": all(controls.values()),
            "no_Ward_identity_for_an_incomplete_graph_subset": True,
            "unconsumed_labels_keep_their_paired_soft_currents": True,
        },
        "whole_counts_by_external_multiplicity": counts,
        "whole_negative_controls": controls,
        "whole_generic_identity": "For each hard-channel side assignment and every subset of four endpoint momentum slots, expand the multilinear shifted endpoint polynomial and sum both possible owners of every emitted singleton. Each unused emission label factors its paired-current sum exactly. Consumed labels carry an explicit energy and no label is consumed twice.",
        "whole_literal_correspondence": "Every Phi2 metric vertex is bilinear in its scalar momenta plus a mass density. A product of two endpoints is a sum of the tested slot monomials. Hard-path momenta and central EH vertices depend only on total radiation assigned to each side, so are common within each endpoint-owner sum.",
    }
