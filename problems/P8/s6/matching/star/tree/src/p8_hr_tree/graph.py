"""Finite oriented-tree incidence and instantaneous nonzero-P components."""

from p8_star.exact import rational


def validate_tree(vertex_count, edges):
    if type(vertex_count) is not int or vertex_count < 1:
        raise ValueError("A finite nonempty integer vertex set is required")
    if not isinstance(edges, (tuple, list)):
        raise TypeError("The finite edge list must be explicit")
    parsed, seen = [], set()
    parent = list(range(vertex_count))
    def representative(vertex):
        while parent[vertex] != vertex:
            vertex = parent[vertex]
        return vertex
    for pair in edges:
        if not isinstance(pair, (tuple, list)) or len(pair) != 2:
            raise ValueError("Each oriented edge has two endpoints")
        i, j = pair
        if any(type(v) is not int or not 0 <= v < vertex_count for v in pair) or i == j:
            raise ValueError("Finite distinct in-range integer endpoints are required")
        key = tuple(sorted(pair))
        if key in seen:
            raise ValueError("Duplicate/reversed duplicate interactions must first be combined")
        seen.add(key)
        left, right = representative(i), representative(j)
        if left == right:
            raise ValueError("Cycles are outside the tree-flux theorem")
        parent[left] = right
        parsed.append((i, j))
    if len(parsed) != vertex_count-1:
        raise ValueError("The supplied graph must be connected; treat forest components separately")
    return tuple(parsed)


def leaf_order(vertex_count, edges):
    parsed = validate_tree(vertex_count, edges)
    adjacent = [set() for _ in range(vertex_count)]
    for index, (i, j) in enumerate(parsed):
        adjacent[i].add(index)
        adjacent[j].add(index)
    remaining = set(range(len(parsed)))
    result = []
    while remaining:
        leaf = next(i for i in range(vertex_count) if len(adjacent[i]) == 1)
        index = next(iter(adjacent[leaf]))
        first, second = parsed[index]
        other = second if leaf == first else first
        result.append((leaf, index, other))
        adjacent[leaf].remove(index)
        adjacent[other].remove(index)
        remaining.remove(index)
    return tuple(result)


def solve_fluxes(vertex_count, edges, divergences):
    """Exact incidence solve; edge flux is + at source and - at target.

    A nonzero supplied divergence is an omission control, not separately
    conserved matter. A tree has only the zero solution for zero data.
    """
    parsed = validate_tree(vertex_count, edges)
    if not isinstance(divergences, (tuple, list)) or len(divergences) != vertex_count:
        raise ValueError("One exact divergence is required per vertex")
    data = list(map(rational, divergences))
    if sum(data) != 0:
        raise ValueError("Incidence data must have zero total divergence")
    result = [rational(0)]*len(parsed)
    for leaf, index, other in leaf_order(vertex_count, parsed):
        sign = 1 if parsed[index][0] == leaf else -1
        result[index] = sign*data[leaf]
        data[other] += data[leaf]
        data[leaf] = rational(0)
    if any(data):
        raise ValueError("The exact leaf-stripping residual did not vanish")
    return tuple(result)


def active_component(vertex_count, edges, polynomial_values, root=0):
    parsed = validate_tree(vertex_count, edges)
    if type(root) is not int or not 0 <= root < vertex_count:
        raise ValueError("A valid distinguished physical vertex is required")
    if not isinstance(polynomial_values, (tuple, list)) or len(polynomial_values) != len(parsed):
        raise ValueError("One exact current P value is required per edge")
    values = tuple(map(rational, polynomial_values))
    found, changed = {root}, True
    while changed:
        changed = False
        for (i, j), p in zip(parsed, values, strict=True):
            if p and ((i in found) != (j in found)):
                found.update((i, j))
                changed = True
    return tuple(sorted(found))


def checks():
    trees = [(1, []), (2, [(0, 1)]), (4, [(0, 1), (1, 2), (2, 3)]),
             (4, [(1, 0), (0, 2), (3, 0)]),
             (6, [(0, 1), (2, 1), (2, 3), (4, 2), (4, 5)])]
    fixtures = []
    for count, edges in trees:
        zero = solve_fluxes(count, edges, [0]*count)
        if any(zero) or len(leaf_order(count, edges)) != count-1:
            raise ValueError("Zero-divergence tree elimination failed")
        chosen = tuple(rational((-1)**i*(i+1)) for i in range(count-1))
        divergence = [rational(0)]*count
        for (i, j), value in zip(edges, chosen, strict=True):
            divergence[i] += value
            divergence[j] -= value
        if solve_fluxes(count, edges, divergence) != chosen:
            raise ValueError("Oriented exact tree incidence replay failed")
        fixtures.append({"vertex_count": count, "edges": edges, "zero_fluxes": list(map(str, zero)),
                         "nonzero_divergence_control_fluxes": list(map(str, chosen))})
    return fixtures
